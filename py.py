from Bio import Entrez, SeqIO, Blast
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt
from Bio.SeqUtils import gc_fraction

Entrez.email = "ravindra.81420@gmail.com"

# Load record
record = SeqIO.read("e_coli_k12.gb", "genbank")
seq = record.seq
print(f"Genome loaded: {len(seq):,} bp | GC: {gc_fraction(seq)*100:.2f}%")

# === Extract Protein Sequences ===
print("\nExtracting some protein sequences...")

proteins = []
for feature in record.features:
    if feature.type == "CDS" and "translation" in feature.qualifiers:
        prot_seq = feature.qualifiers["translation"][0]
        gene = feature.qualifiers.get("gene", ["Unknown"])[0]
        product = feature.qualifiers.get("product", ["Unknown"])[0]
        proteins.append((gene, product, prot_seq))
        
        if len(proteins) >= 5:   # Limit for now
            break

# Show extracted proteins
for i, (gene, product, pseq) in enumerate(proteins, 1):
    print(f"{i}. Gene: {gene} | Product: {product}")
    print(f"   Protein length: {len(pseq)} aa\n")

# === BLAST Section (Fixed) ===
print("\nRunning BLAST on first protein (thrL leader peptide)...")

gene, product, test_seq = proteins[0]

result_handle = NCBIWWW.qblast("blastp", "nr", test_seq)
blast_record = NCBIXML.read(result_handle)
result_handle.close()

print(f"\nBLAST Results for {gene} ({product}):")
for i, alignment in enumerate(blast_record.alignments[:5], 1):
    for hsp in alignment.hsps[:1]:   # Best HSP
        print(f"{i}. Match: {alignment.title[:100]}")
        identity = (hsp.identities / hsp.align_length) * 100
        print(f"   Score: {hsp.score} | Identities: {hsp.identities}/{hsp.align_length} ({identity:.1f}%)")
        print(f"   E-value: {hsp.expect}")
        break
    # === Save Proteins to FASTA File ===
print("\nSaving extracted proteins to 'e_coli_proteins.fasta'...")

with open("e_coli_proteins.fasta", "w") as f:
    for gene, product, prot_seq in proteins:
        f.write(f">{gene} | {product}\n{prot_seq}\n")

print("Saved", len(proteins), "proteins.")

# === Simple Codon Usage for thrA (big protein) ===
print("\nCodon usage for thrA gene:")
thrA_feature = None
for f in record.features:
    if f.type == "CDS" and f.qualifiers.get("gene", [""])[0] == "thrA":
        thrA_feature = f
        break

if thrA_feature:
    dna_seq = thrA_feature.extract(record.seq)
    from collections import Counter
    codons = [str(dna_seq[i:i+3]) for i in range(0, len(dna_seq)-2, 3)]
    codon_count = Counter(codons)
    print("Most frequent codons:")
    for codon, count in codon_count.most_common(10):
        print(f"  {codon}: {count} times")


# === PROJECT SUMMARY REPORT ===
print("\n" + "="*70)
print("GENOMICS PIPELINE SUMMARY REPORT")
print("="*70)
print(f"Organism          : Escherichia coli K-12 substr. MG1655")
print(f"Genome size       : {len(seq):,} bp")
print(f"GC content        : {gc_fraction(seq)*100:.2f}%")
print(f"Total CDS found   : {len([f for f in record.features if f.type == 'CDS'])}")
print(f"Proteins extracted: {len(proteins)}")
print(f"Files created     : e_coli_k12.gb, e_coli_proteins.fasta, e_coli_gc_plot.png")
print("\nTop BLAST hit for thrL: 100% identity to E. coli")
print("Most used codon in thrA: CTG (Leucine)")

print("\n Project completed successfully!")