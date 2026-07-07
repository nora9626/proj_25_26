from fastapi import FastAPI, Depends, HTTPException, status, BackgroundTasks, APIRouter
from fastapi.middleware.cors import CORSMiddleware
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm
from sqlalchemy.orm import Session
from pydantic import BaseModel
from typing import Optional

import datetime as dt
import difflib
import os
import platform
import random
import secrets
import smtplib
import string
import xml.etree.ElementTree as ET
from email.mime.text import MIMEText
from urllib.parse import quote

import bcrypt
import pymysql
import requests
from Bio import Align
from dotenv import load_dotenv

import auth
import database
import models
import schemas
from auth import create_access_token, get_current_user
from database import get_db
from models import Term, UserRole

# Optional dependency. If ChEMBL is not installed, the app will still start.
try:
    from chembl_webresource_client.new_client import new_client
except Exception:
    new_client = None


# ---------------------------------------------------------------------
# Environment and application setup
# ---------------------------------------------------------------------

def smart_load_env():
    current_os = platform.system().lower()
    target_env = "2.env" if current_os == "windows" else ".env"

    if os.path.exists(target_env):
        load_dotenv(target_env)
        print(f"--- Environment loaded from: {target_env} ---")
    else:
        print(f"--- Warning: {target_env} file not found. Using system environment variables. ---")


smart_load_env()

EMAIL_USER = os.getenv("EMAIL_USER")
EMAIL_PASS = os.getenv("EMAIL_PASS")

app = FastAPI(
    title="F8 Hemophilia A Bioinformatics API",
    description=(
        "Educational bioinformatics API for Hemophilia A / F8 gene, "
        "variant classification, sequence comparison, and literature mining."
    ),
    version="1.1.0"
)

router = APIRouter()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


def get_db_connection():
    """
    Direct pymysql connection used by some older sequence-comparison functions.
    SQLAlchemy session is still used for most endpoints.
    """
    db_password = os.getenv("DB_PASS", "")
    return pymysql.connect(
        host=os.getenv("DB_HOST", "localhost"),
        user=os.getenv("DB_USER", "root"),
        password=db_password,
        database=os.getenv("DB_NAME", "dbdictionary"),
        charset="utf8mb4",
        cursorclass=pymysql.cursors.DictCursor
    )


# ---------------------------------------------------------------------
# Request models
# ---------------------------------------------------------------------

class AccessionRequest(BaseModel):
    accession_id: str


class HemophiliaRequest(BaseModel):
    query: str = "Hemophilia A F8 gene therapy"


class LiteratureRequest(BaseModel):
    retmax: int = 5


class F8VariantRequest(BaseModel):
    hgvs: str
    protein: Optional[str] = None
    mutation_type: Optional[str] = None


class DeleteData(BaseModel):
    email: str
    password: str


class SequenceCompareRequest(BaseModel):
    seq1: str
    seq2: str


class TermSaveRequest(BaseModel):
    fasta_seq: str
    disease_class: str
    confidence_score: float


class DnaFullImportRequest(BaseModel):
    accession_id: str
    term_name: Optional[str] = None
    trans: Optional[str] = None
    defe: Optional[str] = None
    disease_class: Optional[str] = "N/A - Genetic Sequence"
    smiles_code: Optional[str] = "N/A"


class DrugRequest(BaseModel):
    drug_name: str


class TermCreate(BaseModel):
    term: str
    trans: str
    defe: str
    user_id: int
    smiles_code: Optional[str] = "N/A"
    picture: Optional[str] = "yyy.jpg"
    status: Optional[str] = "pending"
    sequence: Optional[str] = None
    fasta_seq: Optional[str] = None
    disease_class: Optional[str] = None
    confidence_score: Optional[float] = None


# ---------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------

def translate_to_arabic(text: str) -> str:
    """
    Placeholder translation function.
    Replace later with Google Translate / DeepL / manual translation if needed.
    """
    return text


def is_valid_dna(sequence):
    """
    Validate DNA sequence while allowing N for uncertain bases.
    Returns tuple: (bool, clean_sequence_or_error_message)
    """
    clean_seq = sequence.replace("\n", "").replace(" ", "").upper()
    allowed_chars = set("ACGTN")

    if set(clean_seq).issubset(allowed_chars):
        return True, clean_seq

    return False, "Invalid sequence characters. Allowed DNA characters: A, C, G, T, N."


def send_email(to_email, code):
    if not EMAIL_USER or not EMAIL_PASS:
        print("Error: Email credentials not found in environment variables.")
        return False

    msg = MIMEText(f"كود التحقق الخاص بك هو: {code}")
    msg["Subject"] = "كود استعادة كلمة المرور"
    msg["From"] = EMAIL_USER
    msg["To"] = to_email

    try:
        server = smtplib.SMTP_SSL("smtp.gmail.com", 465)
        server.login(EMAIL_USER, EMAIL_PASS)
        server.sendmail(EMAIL_USER, to_email, msg.as_string())
        server.quit()
        return True
    except Exception as e:
        print(f"Error sending email: {e}")
        return False


def generate_otp_code(length=6):
    return "".join(secrets.choice(string.digits) for _ in range(length))


def verify_password(plain_password, hashed_password):
    return bcrypt.checkpw(
        plain_password[:72].encode("utf-8"),
        hashed_password.encode("utf-8")
    )


def get_password_hash(password):
    return bcrypt.hashpw(
        password[:72].encode("utf-8"),
        bcrypt.gensalt()
    ).decode("utf-8")


def is_admin_user(user: models.User) -> bool:
    role_value = getattr(user.role, "value", user.role)
    return str(role_value) == "admin"


# ---------------------------------------------------------------------
# F8 / Hemophilia A endpoints: the new project direction
# ---------------------------------------------------------------------

@app.post("/api/import-ncbi-f8/")
def import_ncbi_f8(db: Session = Depends(get_db)):
    """
    F8-specific NCBI import for the Hemophilia A project.

    Inputs:
        None.

    Output:
        JSON with NCBI ID, imported term name, and sequence length.

    Assumptions:
        - Educational/research use only.
        - Uses a controlled query focused on F8 in Homo sapiens.
        - Stores the record in the existing terms table.

    Diagnostic checkpoints:
        - Check response status == success.
        - Check sequence_length > 0.
        - Check MySQL terms table for disease_class = 'Hemophilia A / F8'.
    """
    try:
        query = "F8[Gene Name] AND Homo sapiens[Organism]"
        encoded_query = quote(query)

        search_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=nuccore&term={encoded_query}&retmax=5&retmode=json"
        )

        search_res = requests.get(search_url, timeout=10).json()
        id_list = search_res.get("esearchresult", {}).get("idlist", [])

        if not id_list:
            return {
                "status": "error",
                "message": "No F8 sequence records found from NCBI."
            }

        selected_id = id_list[0]

        fasta_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=nuccore&id={selected_id}&rettype=fasta&retmode=text"
        )

        fasta_res = requests.get(fasta_url, timeout=10).text
        lines = fasta_res.splitlines()

        header = lines[0].replace(">", "").strip() if lines else "F8 sequence record"
        clean_sequence = "".join(
            line.strip()
            for line in lines
            if not line.startswith(">")
        )

        new_term = models.Term(
            term=header[:255],
            trans="F8 / Factor VIII sequence record",
            defe=f"NCBI nuccore ID: {selected_id}; controlled query: {query}",
            fasta_seq=clean_sequence,
            sequence=clean_sequence,
            disease_class="Hemophilia A / F8",
            status="approved",
            user_id=46,
            smiles_code="N/A",
            picture="pic/ncbi_logo.png"
        )

        db.add(new_term)
        db.commit()
        db.refresh(new_term)

        return {
            "status": "success",
            "query": query,
            "imported_id": selected_id,
            "term": header,
            "sequence_length": len(clean_sequence)
        }

    except Exception as e:
        db.rollback()
        return {
            "status": "error",
            "message": str(e)
        }


@app.post("/api/import-ncbi/")
def import_ncbi_via_api(db: Session = Depends(get_db)):
    """
    Backward-compatible endpoint.

    The old version searched broadly for 'Hemophilia' and randomly selected one record.
    For the new focused project, this route now redirects to the F8-specific import logic.
    """
    return import_ncbi_f8(db=db)


@app.post("/api/import-hemophilia-f8-literature/")
def import_hemophilia_f8_literature(
    req: LiteratureRequest = LiteratureRequest(),
    db: Session = Depends(get_db)
):
    """
    Retrieves PubMed papers using a controlled Hemophilia A + F8 query.

    Inputs:
        retmax: number of PubMed records to import, capped between 1 and 10.

    Output:
        Imported PMID list and count.

    Filtering:
        The query requires Hemophilia A plus F8/factor VIII plus mutation/variant/gene therapy terms.

    Saved fields:
        term = PMID
        trans = title
        defe = query + abstract
        disease_class = Hemophilia A / F8 literature
    """
    retmax = max(1, min(req.retmax, 10))

    query = 'Hemophilia A AND (F8 OR "factor VIII") AND (mutation OR variant OR "gene therapy")'
    encoded_query = quote(query)

    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=pubmed&term={encoded_query}&retmax={retmax}&retmode=json"
    )

    try:
        search_res = requests.get(search_url, timeout=10).json()
        id_list = search_res.get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"NCBI connection error: {str(e)}"
        )

    imported = []

    for pmid in id_list:
        term_name = f"PMID: {pmid}"

        if db.query(Term).filter(Term.term == term_name).first():
            continue

        fetch_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=pubmed&id={pmid}&retmode=xml"
        )

        try:
            xml_response = requests.get(fetch_url, timeout=10)

            if xml_response.status_code != 200:
                continue

            root = ET.fromstring(xml_response.text)

            title_node = root.find(".//ArticleTitle")
            abstract_nodes = root.findall(".//AbstractText")

            title = "".join(title_node.itertext()) if title_node is not None else "No title"
            abstract = " ".join(
                "".join(node.itertext()).strip()
                for node in abstract_nodes
                if "".join(node.itertext()).strip()
            ) or "No abstract extracted"

            new_term = Term(
                term=term_name,
                trans=title,
                defe=f"Controlled query: {query}\n\nAbstract: {abstract[:1500]}",
                picture="pic/ncbi_logo.png",
                status="approved",
                user_id=46,
                smiles_code="N/A",
                fasta_seq="N/A",
                sequence="N/A",
                disease_class="Hemophilia A / F8 literature"
            )

            db.add(new_term)
            imported.append(pmid)

        except Exception:
            continue

    db.commit()

    return {
        "status": "success",
        "query": query,
        "imported_count": len(imported),
        "pmids": imported
    }


@app.post("/api/import-hemophilia/")
def import_hemophilia_api(
    req: LiteratureRequest = LiteratureRequest(),
    db: Session = Depends(get_db)
):
    """
    Backward-compatible endpoint.

    The old version searched PubMed for a broad 'Hemophilia Gene Therapy' query.
    For the new project, this route now uses the controlled Hemophilia A / F8 literature workflow.
    """
    return import_hemophilia_f8_literature(req=req, db=db)


@app.post("/classify-f8-variant/")
def classify_f8_variant(req: F8VariantRequest):
    """
    Rule-based educational classifier for F8/Hemophilia A variants.

    Inputs:
        hgvs: cDNA or text description, e.g. 'c.1063C>T' or 'intron 22 inversion'
        protein: optional protein change, e.g. 'p.Arg355Ter'

    Output:
        mutation category, interpretation, therapy category, and disclaimer.

    Assumptions:
        - This is not clinical diagnosis.
        - Exact HGVS/ClinVar/EAHAD evidence must be verified before final report.
    """
    hgvs = (req.hgvs or "").strip()
    protein = (req.protein or "").strip()
    combined = f"{hgvs} {protein}".lower()

    if "intron 22" in combined or "intron 1" in combined or "inversion" in combined:
        variant_type = "Structural inversion"
        interpretation = (
            "Large F8 structural variant. Intron 22 and intron 1 inversions "
            "are important causes of severe Hemophilia A and require inversion-specific testing."
        )
        therapy_category = (
            "AAV-mediated gene addition as broad category; gene editing only as future/research concept."
        )

    elif "+" in hgvs or "splice" in combined:
        variant_type = "Splice-site variant"
        interpretation = (
            "Potential disruption of mRNA splicing. Exact effect requires transcript-level evidence."
        )
        therapy_category = (
            "Gene addition as broad category; ASO/splicing strategies only as future/research concept."
        )

    elif "del" in combined or "fs" in combined or "frameshift" in combined or "ins" in combined:
        variant_type = "Frameshift / small insertion-deletion"
        interpretation = (
            "Likely frameshift and premature truncation depending on the variant position."
        )
        therapy_category = (
            "Gene addition as broad category; exon skipping only if specific exon-level evidence exists."
        )

    elif "ter" in combined or "*" in combined or "stop" in combined or "nonsense" in combined:
        variant_type = "Nonsense variant"
        interpretation = (
            "Premature stop codon; predicted truncated FVIII protein."
        )
        therapy_category = (
            "Gene addition as broad category; readthrough/editing only as research/future concept."
        )

    elif ">" in hgvs or "missense" in combined:
        variant_type = "Missense variant"
        interpretation = (
            "Single nucleotide substitution that may cause an amino-acid change. "
            "Functional effect requires evidence or prediction tools."
        )
        therapy_category = (
            "Gene addition as broad category; structural/stability prediction as future work."
        )

    else:
        variant_type = "Unclassified educational entry"
        interpretation = (
            "Insufficient information. Verify HGVS name and source database record."
        )
        therapy_category = "No therapy category assigned until verified."

    return {
        "gene": "F8",
        "disease": "Hemophilia A",
        "hgvs": hgvs,
        "protein": protein,
        "variant_type": variant_type,
        "interpretation": interpretation,
        "therapy_category": therapy_category,
        "disclaimer": "Educational only; not clinical diagnosis or treatment recommendation."
    }


# ---------------------------------------------------------------------
# NCBI gene sequence analysis and suggestions
# ---------------------------------------------------------------------

@app.post("/analyze-gene/")
def analyze_gene(payload: dict):
    """
    Retrieve a nucleotide sequence from NCBI using a gene name, then validate DNA characters.
    Use F8 as the recommended default for this project.
    """
    gene_name = payload.get("gene_name", "F8").strip()

    if not gene_name:
        raise HTTPException(status_code=400, detail="يجب إدخال اسم الجين، مثل F8.")

    query = f"{gene_name}[Gene Name] AND Homo sapiens[Organism]"
    encoded_query = quote(query)

    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=nuccore&term={encoded_query}&retmax=1&retmode=json"
    )

    search_res = requests.get(search_url, timeout=10).json()
    id_list = search_res.get("esearchresult", {}).get("idlist", [])

    if not id_list:
        alt_search_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=nuccore&term={quote(gene_name)}&retmax=1&retmode=json"
        )
        alt_search_res = requests.get(alt_search_url, timeout=10).json()
        id_list = alt_search_res.get("esearchresult", {}).get("idlist", [])

    if not id_list:
        raise HTTPException(
            status_code=404,
            detail="عذراً، لم يتم العثور على تسلسل جيني مطابق لهذا الاسم في NCBI."
        )

    accession = id_list[0]

    fasta_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    )

    fasta_response = requests.get(fasta_url, timeout=10)

    if fasta_response.status_code != 200 or not fasta_response.text:
        raise HTTPException(status_code=500, detail="تعذر جلب التسلسل من قاعدة بيانات NCBI.")

    fasta_text = fasta_response.text
    lines = fasta_text.splitlines()
    raw_sequence = "".join([line.strip() for line in lines if not line.startswith(">")])

    is_valid, clean_seq = is_valid_dna(raw_sequence)

    if not is_valid:
        raise HTTPException(status_code=400, detail="التسلسل المسترجع يحتوي على رموز غير صالحة.")

    return {
        "status": "success",
        "gene_name": gene_name,
        "accession": accession,
        "clean_sequence": clean_seq,
        "sequence_preview": clean_seq[:200],
        "sequence_length": len(clean_seq),
        "disease_classification": "Hemophilia A / F8" if gene_name.upper() == "F8" else "Genetic Research"
    }


@app.get("/suggest-genes/")
def suggest_genes(term: str):
    if len(term) < 2:
        return []

    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=gene&term={quote(term)}*&retmax=5&retmode=json"
    )

    search_res = requests.get(search_url, timeout=10).json()
    id_list = search_res.get("esearchresult", {}).get("idlist", [])

    if not id_list:
        return []

    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=gene&id={','.join(id_list)}&retmode=json"
    )

    summary_res = requests.get(summary_url, timeout=10).json()
    suggestions = []

    for gid in id_list:
        gene_info = summary_res.get("result", {}).get(str(gid))

        if gene_info:
            suggestions.append({
                "id": gid,
                "symbol": gene_info.get("name"),
                "description": gene_info.get("description")
            })

    return suggestions


# ---------------------------------------------------------------------
# Sequence validation and comparison
# ---------------------------------------------------------------------

Nucleotides = ["A", "G", "T", "C"]


def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()

    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False

    return tmpseq


def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}

    for nuc in seq:
        if nuc in tmpFreqDict:
            tmpFreqDict[nuc] += 1

    return tmpFreqDict


@app.post("/analyze-sequence/")
def analyze_sequence(payload: dict):
    user_seq = payload.get("sequence", "").strip().upper()

    if not user_seq:
        raise HTTPException(status_code=400, detail="لم يتم إدخال أي تسلسل.")

    valid_seq = validateSeq(user_seq)

    if not valid_seq:
        raise HTTPException(
            status_code=400,
            detail="التسلسل غير صالح. يجب أن يحتوي فقط على A, C, T, G."
        )

    frequencies = countNucFrequency(valid_seq)

    return {
        "status": "success",
        "validated_sequence": valid_seq,
        "nuc_frequencies": frequencies,
        "message": "تم فحص التسلسل المتغير بنجاح."
    }


@app.post("/compare-two-sequences/")
def compare_two_sequences(payload: dict):
    seq1 = payload.get("seq1", "").replace("\n", "").replace(" ", "").upper()
    seq2 = payload.get("seq2", "").replace("\n", "").replace(" ", "").upper()

    if not seq1 or not seq2:
        raise HTTPException(status_code=400, detail="الرجاء إدخال التسلسلين للمقارنة.")

    aligner = Align.PairwiseAligner()
    alignments = aligner.align(seq1, seq2)

    if not alignments:
        return {
            "status": "success",
            "alignment_score": 0,
            "sequence_1_aligned": seq1,
            "sequence_2_aligned": seq2,
            "message": "لا يوجد تطابق عالي."
        }

    best_alignment = alignments[0]
    score = best_alignment.score

    alignment_str = str(best_alignment).split("\n")
    aligned_1 = alignment_str[0] if len(alignment_str) > 0 else seq1
    aligned_2 = alignment_str[2] if len(alignment_str) > 2 else seq2

    return {
        "status": "success",
        "alignment_score": score,
        "sequence_1_aligned": aligned_1,
        "sequence_2_aligned": aligned_2
    }


@app.post("/generate-random-seq/")
def generate_random_seq(payload: dict):
    length = int(payload.get("length", 3000))
    rndDNAStr = "".join([random.choice(Nucleotides) for _ in range(length)])
    return analyze_sequence({"sequence": rndDNAStr})


@app.get("/generate-random-dna/")
def generate_random_dna():
    rnd_seq = "".join([random.choice(Nucleotides) for _ in range(20)])
    return {"status": "success", "sequence": rnd_seq}


@app.post("/find-closest/")
def find_closest(payload: dict):
    input_seq = payload.get("input_seq", "").strip().upper()
    seq1 = payload.get("seq1", "").strip().upper()
    seq2 = payload.get("seq2", "").strip().upper()

    aligner = Align.PairwiseAligner()

    score1 = 0
    alignments1 = aligner.align(input_seq, seq1)

    if alignments1:
        score1 = alignments1[0].score

    score2 = 0
    alignments2 = aligner.align(input_seq, seq2)

    if alignments2:
        score2 = alignments2[0].score

    if score1 == 0 and score2 == 0:
        winner = "لا يوجد تطابق مع أي من السلسلتين."
    elif score1 >= score2:
        winner = "السلسلة الأولى (Seq 1) هي الأقرب."
    else:
        winner = "السلسلة الثانية (Seq 2) هي الأقرب."

    return {
        "status": "success",
        "score_seq1": score1,
        "score_seq2": score2,
        "closest": winner
    }


@app.post("/compare-with-ncbi/")
def compare_with_ncbi(payload: dict):
    ncbi_seq = payload.get("ncbi_seq", "").strip().upper()
    seq1 = payload.get("seq1", "").strip().upper()
    seq2 = payload.get("seq2", "").strip().upper()

    aligner = Align.PairwiseAligner()

    score1 = 0
    alignments1 = aligner.align(ncbi_seq, seq1)

    if alignments1:
        score1 = alignments1[0].score

    score2 = 0
    alignments2 = aligner.align(ncbi_seq, seq2)

    if alignments2:
        score2 = alignments2[0].score

    if score1 == 0 and score2 == 0:
        winner = "لا يوجد تطابق مع تسلسل NCBI."
    elif score1 >= score2:
        winner = "السلسلة الأولى (Seq 1) هي الأقرب لتسلسل NCBI."
    else:
        winner = "السلسلة الثانية (Seq 2) هي الأقرب لتسلسل NCBI."

    return {
        "status": "success",
        "score_seq1": score1,
        "score_seq2": score2,
        "closest": winner
    }


@app.post("/analyze-sequences/")
def analyze_sequences(req: SequenceCompareRequest):
    s1 = req.seq1.strip().upper()
    s2 = req.seq2.strip().upper()

    def calculate_confidence(seq1, seq2):
        if not seq1 or not seq2:
            return 0.0

        matcher = difflib.SequenceMatcher(None, seq1, seq2)
        return round(matcher.ratio() * 100, 2)

    try:
        connection = get_db_connection()

        with connection.cursor() as cursor:
            query = (
                "SELECT fasta_seq, disease_class "
                "FROM terms "
                "WHERE fasta_seq IS NOT NULL AND fasta_seq != ''"
            )
            cursor.execute(query)
            db_terms = cursor.fetchall()

        connection.close()

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")

    best_match_s1 = {"disease": "لا يوجد تطابق قريب", "score": 0.0}

    if db_terms:
        for term in db_terms:
            fasta = (term.get("fasta_seq") or "").upper()
            score = calculate_confidence(s1, fasta)

            if score > best_match_s1["score"]:
                best_match_s1["score"] = score
                best_match_s1["disease"] = term.get("disease_class") or "Unknown"

    best_match_s2 = {"disease": "لا يوجد تطابق قريب", "score": 0.0}

    if db_terms:
        for term in db_terms:
            fasta = (term.get("fasta_seq") or "").upper()
            score = calculate_confidence(s2, fasta)

            if score > best_match_s2["score"]:
                best_match_s2["score"] = score
                best_match_s2["disease"] = term.get("disease_class") or "Unknown"

    return {
        "status": "success",
        "seq1_analysis": {
            "sequence": s1,
            "associated_disease": best_match_s1["disease"],
            "confidence_score": best_match_s1["score"]
        },
        "seq2_analysis": {
            "sequence": s2,
            "associated_disease": best_match_s2["disease"],
            "confidence_score": best_match_s2["score"]
        }
    }


# ---------------------------------------------------------------------
# Term endpoints
# ---------------------------------------------------------------------

@app.post("/terms/")
def add_term(
    term: schemas.TermSchema,
    db: Session = Depends(database.get_db),
    current_user: models.User = Depends(auth.get_current_user)
):
    term_data = term.dict()

    term_data.setdefault("smiles_code", "N/A")
    term_data.setdefault("fasta_seq", "N/A")
    term_data.setdefault("disease_class", "Genetic Research")

    term_data["user_id"] = current_user.id
    term_data["status"] = "pending"

    new_term = models.Term(**term_data)

    try:
        db.add(new_term)
        db.commit()
        db.refresh(new_term)

    except Exception as e:
        db.rollback()
        print("Database Error:", e)
        raise HTTPException(status_code=500, detail=str(e))

    return {
        "message": "تمت الإضافة بنجاح",
        "term_id": new_term.id
    }


@app.put("/terms/update")
def update_term(data: dict, db: Session = Depends(get_db)):
    term_id = data.get("id")

    if not term_id:
        raise HTTPException(status_code=400, detail="ID مفقود")

    term = db.query(models.Term).filter(models.Term.id == term_id).first()

    if not term:
        raise HTTPException(status_code=404, detail="المصطلح غير موجود")

    term.term = data.get("term", term.term)
    term.trans = data.get("trans", term.trans)
    term.defe = data.get("defe", term.defe)
    term.smiles_code = data.get("smiles_code", term.smiles_code)
    term.fasta_seq = data.get("fasta_seq", term.fasta_seq)
    term.sequence = data.get("sequence", term.sequence)
    term.disease_class = data.get("disease_class", term.disease_class)
    term.confidence_score = data.get("confidence_score", term.confidence_score)

    if "status" in data:
        term.status = data["status"]

    db.commit()
    db.refresh(term)

    return {"status": "success"}


@app.delete("/terms/{term_id}")
def delete_term(term_id: int, db: Session = Depends(database.get_db)):
    term = db.query(models.Term).filter(models.Term.id == term_id).first()

    if not term:
        raise HTTPException(status_code=404, detail="غير موجود")

    db.delete(term)
    db.commit()

    return {"message": "تم الحذف"}


@app.put("/terms/approve/{term_id}")
def approve_term(
    term_id: int,
    db: Session = Depends(database.get_db),
    current_user: models.User = Depends(auth.get_current_user)
):
    if not is_admin_user(current_user):
        raise HTTPException(status_code=403, detail="لا تملك صلاحية الموافقة على المصطلحات")

    term = db.query(models.Term).filter(models.Term.id == term_id).first()

    if not term:
        raise HTTPException(status_code=404, detail="المصطلح غير موجود")

    term.status = "approved"
    db.commit()

    return {"message": "تمت الموافقة على المصطلح وعرضه للجمهور بنجاح"}


@app.get("/terms/uu/")
def get_terms(term_id: Optional[int] = None, db: Session = Depends(database.get_db)):
    if term_id is not None:
        term = db.query(models.Term).filter(models.Term.id == term_id).first()

        if not term:
            raise HTTPException(status_code=404, detail="المصطلح غير موجود")

        return term

    terms = (
        db.query(models.Term)
        .filter(models.Term.status != "rejected")
        .order_by(models.Term.id.desc())
        .all()
    )

    return terms


@app.get("/terms/public/")
def get_public_terms(db: Session = Depends(database.get_db)):
    return db.query(models.Term).filter(models.Term.status == "approved").all()


@app.get("/terms/user/{user_id}")
def get_user_terms(user_id: int, db: Session = Depends(get_db)):
    terms = (
        db.query(models.Term)
        .filter(models.Term.user_id == user_id)
        .order_by(models.Term.id.desc())
        .all()
    )

    return [
        {
            "term": t.term,
            "status": t.status
        }
        for t in terms
    ]


@app.post("/terms/add-from-bot/")
def add_term_from_bot(data: dict, db: Session = Depends(get_db)):
    new_term = models.Term(
        term=data.get("term"),
        trans=data.get("trans"),
        defe=data.get("defe"),
        smiles_code=data.get("smiles_code", "N/A"),
        user_id=data.get("user_id", 46),
        status=data.get("status", "pending"),
        picture=data.get("picture", "pic/ncbi_logo.png"),
        fasta_seq=data.get("fasta_seq", "N/A"),
        sequence=data.get("sequence", "N/A"),
        disease_class=data.get("disease_class", "Genetic Research")
    )

    db.add(new_term)
    db.commit()
    db.refresh(new_term)

    return {
        "status": "success",
        "term_id": new_term.id
    }


@app.get("/terms/count/")
def get_terms_count(db: Session = Depends(get_db)):
    count = db.query(models.Term).count()
    return {"total": count}


@app.get("/terms/")
def get_all_terms(db: Session = Depends(get_db)):
    terms = (
        db.query(models.Term)
        .filter(models.Term.status != "rejected")
        .order_by(models.Term.id.desc())
        .all()
    )

    return terms


# ---------------------------------------------------------------------
# Authentication and user account endpoints
# ---------------------------------------------------------------------

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token")
otp_memory = {}


@app.post("/login/")
def login(
    form_data: OAuth2PasswordRequestForm = Depends(),
    db: Session = Depends(database.get_db)
):
    user = (
        db.query(models.User)
        .filter(models.User.username == form_data.username)
        .first()
    )

    if not user or not verify_password(form_data.password, user.passwor):
        raise HTTPException(status_code=400, detail="اسم المستخدم أو كلمة المرور غير صحيحة")

    access_token = create_access_token(data={"sub": user.username})

    return {
        "access_token": access_token,
        "token_type": "bearer"
    }


@app.post("/otp/send-code/")
def send_otp_code(
    data: dict,
    background_tasks: BackgroundTasks,
    db: Session = Depends(database.get_db)
):
    email = data.get("email")
    code = str(random.randint(100000, 999999))

    otp_memory[code] = email

    expire_time = dt.datetime.now(dt.timezone.utc) + dt.timedelta(minutes=10)
    new_otp = models.OTP(code=code, expires_at=expire_time)

    db.add(new_otp)
    db.commit()

    background_tasks.add_task(send_email, email, code)

    return {"message": "تم إرسال الكود"}


@app.post("/users/finalize-signup/")
def finalize_signup(data: dict, db: Session = Depends(database.get_db)):
    code = data.get("code")

    otp_record = db.query(models.OTP).filter(models.OTP.code == code).first()
    email = otp_memory.get(code)

    if not otp_record or not email:
        raise HTTPException(status_code=400, detail="❌ الكود غير صحيح أو انتهت صلاحيته.")

    new_user = models.User(
        username=data.get("username"),
        passwor=get_password_hash(data.get("password")),
        email=email,
        phone=data.get("phone")
    )

    db.add(new_user)
    db.delete(otp_record)

    if code in otp_memory:
        del otp_memory[code]

    db.commit()

    return {
        "status": "success",
        "message": "تم إنشاء الحساب بنجاح"
    }


@app.post("/user/delete/")
def delete_user(email: str, code: str, db: Session = Depends(database.get_db)):
    user = db.query(models.User).filter(models.User.email == email).first()

    if not user:
        raise HTTPException(status_code=404, detail="مستخدم غير موجود")

    otp_record = (
        db.query(models.OTP)
        .filter(
            models.OTP.user_id == user.id,
            models.OTP.delete_otp == code
        )
        .first()
    )

    if not otp_record:
        raise HTTPException(status_code=400, detail="الكود غير صحيح")

    db_time = otp_record.delete_expire

    if db_time.tzinfo is None:
        db_time = db_time.replace(tzinfo=dt.timezone.utc)

    if db_time < dt.datetime.now(dt.timezone.utc):
        raise HTTPException(status_code=400, detail="الكود منتهي الصلاحية")

    db.delete(user)
    db.delete(otp_record)
    db.commit()

    return {"status": "deleted"}


@app.post("/otp/send/")
def send_otp(
    email: str,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db)
):
    user = db.query(models.User).filter(models.User.email == email).first()

    if not user:
        raise HTTPException(status_code=404, detail="مستخدم غير موجود")

    otp_code = generate_otp_code(6)
    background_tasks.add_task(send_email, email, otp_code)

    new_otp = models.OTP(
        user_id=user.id,
        reset_code=otp_code,
        expires_at=dt.datetime.now(dt.timezone.utc) + dt.timedelta(minutes=15)
    )

    db.add(new_otp)
    db.commit()

    return {"message": "تم إرسال الكود، قد يستغرق وصوله ثوانٍ معدودة."}


@app.post("/otp/verify/")
def verify_otp(
    email: str,
    code: str,
    db: Session = Depends(database.get_db)
):
    otp_record = (
        db.query(models.OTP)
        .join(models.User)
        .filter(
            models.User.email == email,
            models.OTP.reset_code == code,
            models.OTP.expires_at > dt.datetime.now(dt.timezone.utc)
        )
        .first()
    )

    if not otp_record:
        print(f"DEBUG: No OTP found for {email} with code {code}")
        raise HTTPException(status_code=400, detail="الكود غير صحيح أو منتهي الصلاحية")

    db.delete(otp_record)
    db.commit()

    return {
        "status": "success",
        "message": "تم التحقق بنجاح"
    }


@app.post("/otp/send-delete/")
def send_delete_otp(
    data: DeleteData,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db)
):
    user = db.query(models.User).filter(models.User.email == data.email).first()

    if not user or not verify_password(data.password, user.passwor):
        raise HTTPException(status_code=400, detail="كلمة المرور غير صحيحة")

    otp_code = generate_otp_code(6)
    expire_time = dt.datetime.now(dt.timezone.utc) + dt.timedelta(minutes=5)

    new_otp = models.OTP(
        user_id=user.id,
        delete_otp=otp_code,
        delete_expire=expire_time
    )

    db.add(new_otp)
    db.commit()

    background_tasks.add_task(send_email, data.email, otp_code)

    return {
        "message": "تم إرسال الكود",
        "expires_at": expire_time.isoformat()
    }


@app.get("/users/me/")
def get_me(current_user: models.User = Depends(auth.get_current_user)):
    return {
        "username": current_user.username,
        "email": current_user.email,
        "phone": current_user.phone
    }


@app.get("/users/mee/")
def get_user_profile(current_user: models.User = Depends(get_current_user)):
    role_str = (
        str(current_user.role.value)
        if hasattr(current_user.role, "value")
        else str(current_user.role)
    )

    response_data = {
        "username": current_user.username,
        "email": current_user.email,
        "phone": current_user.phone,
        "role": role_str
    }

    print(f"DEBUG: Data sent to PHP: {response_data}")

    return response_data


@app.post("/user/request-update/")
def request_user_update(
    data: dict,
    background_tasks: BackgroundTasks,
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    new_email = data.get("new_email")
    new_phone = data.get("new_phone")
    new_name = data.get("new_name")

    code_val = generate_otp_code(6)
    expire_time = dt.datetime.now(dt.timezone.utc) + dt.timedelta(minutes=5)

    new_otp = models.OTP(
        user_id=current_user.id,
        edit_code=code_val,
        edit_expire=expire_time
    )

    db.add(new_otp)

    current_user.pending_email = new_email
    current_user.pending_phone = new_phone
    current_user.pending_name = new_name

    db.commit()

    background_tasks.add_task(send_email, new_email, code_val)

    return {"message": "تم إرسال كود التأكيد إلى إيميلك الجديد"}


@app.post("/otp/change/")
def change_otp(
    data: dict,
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    input_code = data.get("code")

    otp_record = (
        db.query(models.OTP)
        .filter(
            models.OTP.user_id == current_user.id,
            models.OTP.edit_code == input_code
        )
        .first()
    )

    if not otp_record:
        raise HTTPException(status_code=400, detail="الكود غير موجود")

    if otp_record.edit_expire < dt.datetime.now(dt.timezone.utc).replace(tzinfo=None):
        db.delete(otp_record)
        db.commit()
        raise HTTPException(status_code=400, detail="الكود منتهي الصلاحية")

    if current_user.pending_email:
        current_user.email = current_user.pending_email

    if current_user.pending_phone:
        current_user.phone = current_user.pending_phone

    if current_user.pending_name:
        current_user.username = current_user.pending_name

    current_user.pending_email = None
    current_user.pending_phone = None
    current_user.pending_name = None

    db.delete(otp_record)
    db.commit()

    return {
        "status": "success",
        "message": "تم التحديث بنجاح"
    }


@app.post("/password/request-reset/")
def request_password_reset(
    email: str,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db)
):
    user = db.query(models.User).filter(models.User.email == email).first()

    if not user:
        raise HTTPException(status_code=404, detail="هذا البريد غير مسجل")

    otp_code = generate_otp_code(6)
    expire_time = dt.datetime.now(dt.timezone.utc) + dt.timedelta(minutes=10)

    new_otp = models.OTP(
        user_id=user.id,
        reset_code=otp_code,
        reset_expire=expire_time
    )

    db.add(new_otp)
    db.commit()

    background_tasks.add_task(send_email, email, otp_code)

    return {"message": "تم إرسال كود إعادة تعيين كلمة المرور"}


@app.post("/password/verify-code/")
def verify_code(data: dict, db: Session = Depends(get_db)):
    email = data.get("email")
    code = data.get("code")

    user = db.query(models.User).filter(models.User.email == email).first()

    if not user:
        raise HTTPException(status_code=404, detail="هذا البريد غير مسجل")

    otp_record = (
        db.query(models.OTP)
        .filter(
            models.OTP.user_id == user.id,
            models.OTP.reset_code == code,
            models.OTP.reset_expire >= dt.datetime.now(dt.timezone.utc)
        )
        .first()
    )

    if not otp_record:
        raise HTTPException(status_code=400, detail="الكود غير صحيح أو منتهي الصلاحية")

    return {
        "status": "success",
        "message": "الكود صحيح"
    }


@app.post("/password/reset/")
def reset_password(data: dict, db: Session = Depends(get_db)):
    email = data.get("email")
    new_password = data.get("new_password")

    user = db.query(models.User).filter(models.User.email == email).first()

    if not user:
        raise HTTPException(status_code=404, detail="هذا البريد غير مسجل")

    user.passwor = get_password_hash(new_password)
    db.commit()

    return {"message": "تم تغيير كلمة المرور بنجاح"}


@app.post("/password/modifypa-otp/")
def modifypa_otp(data: dict, db: Session = Depends(get_db)):
    email = data.get("email")
    otp_input = data.get("code")

    user = db.query(models.User).filter(models.User.email == email).first()

    if not user:
        raise HTTPException(status_code=404, detail="مستخدم غير موجود")

    otp_record = (
        db.query(models.OTP)
        .filter(
            models.OTP.user_id == user.id,
            models.OTP.reset_code == otp_input
        )
        .first()
    )

    if not otp_record:
        raise HTTPException(status_code=400, detail="❌ الكود غير صحيح")

    expire_time = otp_record.reset_expire

    if expire_time.tzinfo is None:
        expire_time = expire_time.replace(tzinfo=dt.timezone.utc)

    if expire_time < dt.datetime.now(dt.timezone.utc):
        raise HTTPException(status_code=400, detail="❌ الكود منتهي الصلاحية")

    return {
        "status": "success",
        "message": "الكود صحيح"
    }


# ---------------------------------------------------------------------
# ChEMBL drug import
# ---------------------------------------------------------------------

@app.post("/api/import-dna-complete/")
def add_drug_via_api(req: DrugRequest, db: Session = Depends(get_db)):
    """
    Legacy ChEMBL import endpoint.
    Kept for compatibility, but it is not central to the F8/Hemophilia A MVP.
    """
    if new_client is None:
        raise HTTPException(
            status_code=500,
            detail="chembl_webresource_client is not installed in this environment."
        )

    drug_name = req.drug_name.strip()

    if not drug_name:
        raise HTTPException(status_code=400, detail="يجب إدخال اسم الدواء.")

    molecule = new_client.molecule
    res = molecule.search(drug_name)

    if not res:
        raise HTTPException(
            status_code=404,
            detail=f"لم يتم العثور على الدواء: {drug_name} في ChEMBL."
        )

    data = res[0]
    name_en = data.get("pref_name", drug_name)

    structures = data.get("molecule_structures")
    smiles = structures.get("canonical_smiles", "N/A") if structures else "N/A"

    properties = data.get("molecule_properties")
    mw = properties.get("full_mwt", "Unknown") if properties else "Unknown"

    description = f"Molecular Weight: {mw} g/mol.\nSource: ChEMBL Database."

    try:
        db_term = Term(
            term=name_en,
            trans=name_en,
            defe=description,
            smiles_code=smiles,
            picture="pic/chembl_logo.png",
            status="approved",
            user_id=46,
            fasta_seq="N/A",
            sequence="N/A",
            disease_class="Hemophilia Related / ChEMBL"
        )

        db.add(db_term)
        db.commit()
        db.refresh(db_term)

    except Exception as e:
        db.rollback()
        raise HTTPException(
            status_code=500,
            detail=f"خطأ في حفظ الدواء بقاعدة البيانات: {str(e)}"
        )

    return {
        "status": "success",
        "message": f"تم جلب بيانات {name_en} وحفظها بنجاح.",
        "drug_name": name_en,
        "smiles": smiles
    }


# ---------------------------------------------------------------------
# Local run
# ---------------------------------------------------------------------

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
