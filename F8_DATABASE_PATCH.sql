-- F8 Hemophilia A Variant Explorer database patch
-- Purpose: add backend fields expected by the updated FastAPI main.py.
-- Run this once in phpMyAdmin or MySQL before testing the F8 endpoints.
-- This patch does not delete, rename, or move any existing data.

ALTER TABLE terms
ADD COLUMN sequence LONGTEXT NULL,
ADD COLUMN fasta_seq LONGTEXT NULL,
ADD COLUMN disease_class VARCHAR(255) NULL,
ADD COLUMN confidence_score DECIMAL(5,2) NULL,
ADD COLUMN add_date DATETIME DEFAULT CURRENT_TIMESTAMP;
-- Allow longer NCBI FASTA record names / headers
ALTER TABLE terms MODIFY COLUMN term VARCHAR(255) NOT NULL;

-- Diagnostic checkpoint: run after the ALTER TABLE command
-- DESCRIBE terms;

-- Optional diagnostic query after importing F8 records
-- SELECT id, term, trans, disease_class, LENGTH(fasta_seq) AS fasta_length
-- FROM terms
-- WHERE disease_class LIKE '%F8%' OR term LIKE '%F8%' OR trans LIKE '%F8%'
-- ORDER BY id DESC
-- LIMIT 10;
