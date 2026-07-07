# F8 Hemophilia A Variant Explorer - Testing Steps

## What changed

The backend is refocused on an educational Hemophilia A / F8 workflow instead of broad random Hemophilia retrieval.

Main endpoints added or refocused:

- `POST /api/import-ncbi-f8/`
  - Retrieves an F8-related sequence record using a controlled NCBI query.
  - Query: `F8[Gene Name] AND Homo sapiens[Organism]`

- `POST /api/import-ncbi/`
  - Kept for compatibility, but redirected to the safer F8-specific import workflow.

- `POST /api/import-hemophilia-f8-literature/`
  - Retrieves PubMed records using a controlled Hemophilia A / F8 query.

- `POST /classify-f8-variant/`
  - Educational rule-based classifier for F8 variant descriptions.
  - It does not diagnose and does not recommend treatment.

## Before testing

Run the database patch once in MySQL/phpMyAdmin:

```sql
ALTER TABLE terms
ADD COLUMN sequence LONGTEXT NULL,
ADD COLUMN fasta_seq LONGTEXT NULL,
ADD COLUMN disease_class VARCHAR(255) NULL,
ADD COLUMN confidence_score DECIMAL(5,2) NULL,
ADD COLUMN add_date DATETIME DEFAULT CURRENT_TIMESTAMP;
```

Diagnostic check:

```sql
DESCRIBE terms;
```

Expected fields should include:

- `sequence`
- `fasta_seq`
- `disease_class`
- `confidence_score`
- `add_date`

## Run backend

From the project folder:

```bash
uvicorn main:app --reload
```

Open:

```text
http://127.0.0.1:8000/docs
```

## Test 1: F8 NCBI import

Endpoint:

```text
POST /api/import-ncbi-f8/
```

Expected output contains:

```json
{
  "status": "success",
  "query": "F8[Gene Name] AND Homo sapiens[Organism]",
  "sequence_length": 12345
}
```

The exact sequence length may differ depending on the NCBI record returned.

Database check:

```sql
SELECT id, term, trans, disease_class, LENGTH(fasta_seq) AS fasta_length
FROM terms
ORDER BY id DESC
LIMIT 5;
```

## Test 2: Hemophilia A / F8 PubMed literature import

Endpoint:

```text
POST /api/import-hemophilia-f8-literature/
```

Example JSON body:

```json
{
  "retmax": 5
}
```

Expected output:

```json
{
  "status": "success",
  "imported_count": 1,
  "pmids": ["..."]
}
```

The imported count can be 0 if the papers already exist in the `terms` table.

## Test 3: Variant classifier

Endpoint:

```text
POST /classify-f8-variant/
```

Example 1:

```json
{
  "hgvs": "c.1063C>T",
  "protein": "p.Arg355Ter"
}
```

Expected classification:

```text
Nonsense variant
```

Example 2:

```json
{
  "hgvs": "intron 22 inversion",
  "protein": ""
}
```

Expected classification:

```text
Structural inversion
```

Example 3:

```json
{
  "hgvs": "c.143+1G>A",
  "protein": ""
}
```

Expected classification:

```text
Splice-site variant
```

## What to say in the presentation

Use this wording:

> The implemented system uses controlled F8/Hemophilia A queries to retrieve records from NCBI and PubMed. Retrieved records are parsed, stored in the website database, and displayed to the user. A simple educational rule-based classifier maps F8 variant descriptions to broad mutation categories such as nonsense, missense, splice-site, frameshift, or inversion. The system is for educational/research support only and does not provide diagnosis or treatment recommendations.

Avoid saying:

> The system develops genetic drugs.

Better wording:

> The system links mutation categories to research-level therapeutic strategy concepts.
