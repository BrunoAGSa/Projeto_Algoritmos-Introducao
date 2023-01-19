select * from Sequences;
select * from Locus;
select * from Genbank;
select * from Features;
select * from CDS;
select * from Authors;
select * from reference_authors;
select * from Pubmed_Info;
select * from genbank_reference;
select * from reference;

# Exemplos de uso 
# 1 Quais os artigos Pubmed (DOI) relacioanados ao Organismo Ectocarpus siliculosus

select DOI from Pubmed_Info
join reference on Pubmed_Info.ID_journal = reference.Journal_ID
join genbank_reference on reference.Journal_ID = genbank_reference.journal_connect
join Genbank on genbank_reference.Id_version_connect = Genbank.ID_version_seq
where Organism = "Ectocarpus siliculosus";

# 2 Foi analisada uma Proteina cujo nome é DnaA-like replication initiation protein.
# Qual o Organismo a qual esta pertence?

select Organism from Genbank
join Features on Genbank.ID_version_seq = Features.ID_version_genbank
join CDS on Features.ID_version_genbank = CDS.Id_version_features
where Protein = '[''DnaA-like replication initiation protein'']'

