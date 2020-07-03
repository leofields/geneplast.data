#!/bin/bash --
gzip -d odb10v1_*.tab.gz
sqlite3 -batch $1 <<EOF
CREATE TABLE odb_species (tax_id INTEGER, org_id TEXT, scientific_name TEXT, genome_assembly_id TEXT, count_genes INTEGER, count_ogs INTEGER, mapping_type TEXT);
CREATE TABLE odb_orthgroups (og_id TEXT, level_tax_id INTEGER, og_name TEXT);
CREATE TABLE odb_genes (gene_id TEXT, org_id TEXT, protein_id TEXT, synonyms_1 TEXT, synonyms_2 TEXT, synonyms_3 INTEGER, description TEXT);
CREATE TABLE odb_orthgroups_genes (og_id TEXT, gene_id TEXT);
.separator "\t"
.import odb10v1_species.tab odb_species
.import odb10v1_OGs.tab odb_orthgroups
.import odb10v1_genes.tab odb_genes
.import odb10v1_OG2genes.tab odb_orthgroups_genes
.headers on
.mode csv
.once orthodb_cogdata.csv
SELECT protein_id, tax_id, og_id FROM odb_orthgroups JOIN odb_orthgroups_genes USING (og_id) JOIN odb_genes USING (gene_id) JOIN odb_species USING (org_id) WHERE level_tax_id = 2759;
EOF
