# Biopython-pipeline-for-protein-and-nucleotide
This pipeline explores the various biological functions peroformed on nucleotide and protein using biopython. 
# Biopython-assignment (enhanced)
# Prapti Soni
# NOTE: run this in an environment with internet access and Biopython installed.

 
!pip install --upgrade biopython

from Bio import Entrez, SeqIO, SwissProt, Medline, Phylo
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, seq1
from Bio.PDB import PDBList, PDBParser, PDBIO, MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder, is_aa
from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import motifs
import xml.etree.ElementTree as ET
import os

# -------------------------------
# Settings
# -------------------------------
Entrez.email = "praptisoni204@gmail.com"  # CHANGE to your email
NUC_ACC = "NM_001355006"
PDB_ID = "1TUP"
OUTDIR = "pipeline_output"
os.makedirs(OUTDIR, exist_ok=True)

# -------------------------------
# Fetch nucleotide sequence (GenBank + FASTA parsing)
# -------------------------------
print("\n--- Fetching nucleotide (GenBank + FASTA) ---")
try:
    # GenBank record
    handle_gb = Entrez.efetch(db="nucleotide", id=NUC_ACC, rettype="gb", retmode="text")
    gb_record = SeqIO.read(handle_gb, "genbank")
    handle_gb.close()
    print("GenBank ID:", gb_record.id)
    print("Description:", gb_record.description)
    # Save GenBank locally
    gb_path = os.path.join(OUTDIR, f"{NUC_ACC}.gb")
    SeqIO.write(gb_record, gb_path, "genbank")

    # FASTA (same accession)
    handle_fa = Entrez.efetch(db="nucleotide", id=NUC_ACC, rettype="fasta", retmode="text")
    fa_record = SeqIO.read(handle_fa, "fasta")
    handle_fa.close()
    nuc_seq = fa_record.seq
except Exception as e:
    print("Failed to fetch nucleotide record:", e)
    nuc_seq = Seq("")

if nuc_seq:
    protein_from_nuc = nuc_seq.translate(to_stop=True)
    print("Length of nucleotide:", len(nuc_seq))
    print("GC%: {:.2f}%".format(gc_fraction(nuc_seq) * 100))
else:
    protein_from_nuc = Seq("")

# -------------------------------
# Fetch Swiss-Prot (UniProt/SwissProt) sequence
# -------------------------------
# This example uses Bio.SwissProt to parse a local/uniprot-formatted record or ExPASy if available.
print("\n--- Swiss-Prot fetching/parsing example ---")
try:
    # Example: fetch raw swissprot entry via Entrez is not appropriate; normally use ExPASy.get_sprot_raw
    # We'll attempt to parse a local UniProt file if present, otherwise demo placeholder
    swp_path = os.path.join(OUTDIR, "uniprot_example.dat")
    if os.path.exists(swp_path):
        with open(swp_path) as handle:
            for record in SwissProt.parse(handle):
                print("SwissProt accession:", record.accessions)
                break
    else:
        print("No local SwissProt file found at", swp_path)
        print("If you want to fetch SwissProt/UniProt programmatically, use Bio.ExPASy and Bio.SwissProt.parse on the downloaded file.")
except Exception as e:
    print("SwissProt parsing failed:", e)

# -------------------------------
# Fetch protein structure (PDB/mmCIF) and read/write
# -------------------------------
print("\n--- PDB/mmCIF fetch and read/write ---")
try:
    pdbl = PDBList()
    # retrieve_pdb_file returns a path to the downloaded file (PDB format by default)
    pdb_file = pdbl.retrieve_pdb_file(PDB_ID, pdir=OUTDIR, file_format="pdb")
    print("Downloaded PDB to:", pdb_file)

    # Parse using PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, pdb_file)

    # Extract sequence from first chain similar to previous code
    ppb = PPBuilder()
    prot_seq = ""
    for model in structure:
        for chain in model:
            seq_list = []
            for pp in ppb.build_peptides(chain):
                for res in pp:
                    if is_aa(res, standard=True):
                        seq_list.append(seq1(res.get_resname()))
            if seq_list:
                prot_seq = "".join(seq_list)
                print(f"Chain {chain.id} length: {len(prot_seq)}")
                break
        if prot_seq:
            break

    # Write structure back out (example: write as PDB)
    io = PDBIO()
    io.set_structure(structure)
    out_pdb_path = os.path.join(OUTDIR, f"{PDB_ID}_copy.pdb")
    io.save(out_pdb_path)
    print("Saved copy of structure to:", out_pdb_path)

    # If mmCIF is desired and available, show how to parse mmCIF
    cif_path = pdbl.retrieve_pdb_file(PDB_ID, pdir=OUTDIR, file_format="mmCif")
    if os.path.exists(cif_path):
        mmcif_parser = MMCIFParser()
        mmcif_struct = mmcif_parser.get_structure(PDB_ID + "_cif", cif_path)
        print("Parsed mmCIF structure (example)")
except Exception as e:
    print("PDB handling failed:", e)

# -------------------------------
# Pairwise Alignment (protein_short vs prot_short)
# -------------------------------
print("\n--- Pairwise alignment ---")
protein_short = protein_from_nuc[:200]
prot_short = prot_seq[:200] if prot_seq else ""
if len(protein_short) == 0 or len(prot_short) == 0:
    print("One of the sequences is empty. Alignment skipped.")
else:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    global_alignment = aligner.align(protein_short, prot_short)
    if len(global_alignment) > 0:
        print(global_alignment[0])
    aligner.mode = "local"
    local_alignment = aligner.align(protein_short, prot_short)
    if len(local_alignment) > 0:
        print(local_alignment[0])

# -------------------------------
# Find first start codon
# -------------------------------
print("\n--- Start codon search ---")
for i in range(0, len(nuc_seq)-2, 3):
    if nuc_seq[i:i+3] == "ATG":
        print("First start codon found at nucleotide position:", i)
        break

# -------------------------------
# BLAST (protein)
# -------------------------------
print("\n--- BLAST example (blastp) ---")
try:
    blast_file = os.path.join(OUTDIR, "blast_results.xml")
    if protein_from_nuc:
        # NCBIWWW.qblast can be slow and requires internet; handle gracefully
        print("Submitting BLAST (blastp) request... this may take several minutes")
        result_handle = NCBIWWW.qblast("blastp", "nr", str(protein_from_nuc), hitlist_size=5)
        with open(blast_file, "w") as f:
            f.write(result_handle.read())
        result_handle.close()

        with open(blast_file) as f:
            blast_records = NCBIXML.read(f)
        print("Top BLAST hits:")
        for i, alignment in enumerate(blast_records.alignments[:5]):
            print(i+1, alignment.hit_id, alignment.length, alignment.hsps[0].expect)
    else:
        print("Protein sequence empty; BLAST skipped.")
except Exception as e:
    print("BLAST failed/was skipped:", e)

# -------------------------------
# Parse MEDLINE (PubMed) records for the GenBank record (example)
# -------------------------------
print("\n--- MEDLINE / PubMed fetch & parse ---")
try:
    # Example: use a PubMed ID if available in the GenBank features (look for PUBMED in references)
    pubmed_ids = []
    if hasattr(gb_record, 'annotations') and 'references' in gb_record.annotations:
        for ref in gb_record.annotations.get('references', []):
            # reference.pubmed_id may or may not be present depending on parsing
            if hasattr(ref, 'pubmed_id') and ref.pubmed_id:
                pubmed_ids.append(ref.pubmed_id)
    # As a fallback, provide an example PubMed ID list
    if not pubmed_ids:
        # Example PubMed ID (replace with real IDs you want)
        pubmed_ids = ["31452104"]  # placeholder

    handle_pm = Entrez.efetch(db="pubmed", id=','.join(pubmed_ids), rettype="medline", retmode="text")
    records = Medline.parse(handle_pm)
    for rec in records:
        print("PMID:", rec.get('PMID'))
        print("Title:", rec.get('TI'))
        break
    handle_pm.close()
except Exception as e:
    print("MEDLINE fetch/parse failed:", e)

# -------------------------------
# GEO / GEORecords example (fetch XML and parse minimal fields)
# -------------------------------
print("\n--- GEO/GDS example (fetch via Entrez and parse XML) ---")
try:
    # Example: fetch a GDS (Geo DataSet) record by GDS id; replace with a real GDS id if desired
    gds_id = "GDS507"  # placeholder example
    handle_gds = Entrez.efetch(db="gds", id=gds_id, rettype="full", retmode="xml")
    xml_text = handle_gds.read()
    handle_gds.close()
    root = ET.fromstring(xml_text)
    # Print some tags present (this depends on the returned XML shape)
    print("GDS XML root tag:", root.tag)
except Exception as e:
    print("GEO/GDS fetch failed or placeholder used:", e)

# -------------------------------
# Phylogenetics: build tree from multiple sequences (Distance-based)
# -------------------------------
print("\n--- Phylogenetic tree from FASTA sequences example ---")
try:
    # For demo create a small FASTA file of related sequences
    fasta_demo = os.path.join(OUTDIR, "demo_seqs.fasta")
    demo_seqs = [
        ("seq1", "MKTAYIAKQRQISFVKSHFSRQDILDLWQ"),
        ("seq2", "MKTAYIAKQRQISFVKSHFSRQDIVDLWQ"),
        ("seq3", "MKTAYIAKQRQISFVKSHFSRQDILDLWK"),
    ]
    with open(fasta_demo, "w") as fh:
        for n, s in demo_seqs:
            fh.write(f">{n}\n{s}\n")

    alignment = AlignIO.read(fasta_demo, "fasta")  # these are equal-length strings in demo
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    tree_path = os.path.join(OUTDIR, "demo_tree.xml")
    Phylo.write(tree, tree_path, "phyloxml")
    print("Wrote phylogenetic tree to:", tree_path)
    print("Tree (ASCII):")
    Phylo.draw_ascii(tree)
except Exception as e:
    print("Phylogenetics example failed:", e)

# -------------------------------
# Sequence motif analysis using Bio.motifs
# -------------------------------
print("\n--- Motif analysis example (Bio.motifs) ---")
try:
    # Using the demo sequences above to create a motif
    seq_records = [Seq(s) for _, s in demo_seqs]
    motif_instances = [m for m in seq_records]
    # Biopython motifs expects Seq objects and a list of sequences
    m = motifs.create(seq_records)
    print("Motif length:", len(m))
    # Calculate consensus
    try:
        print("Consensus:", m.consensus)
    except Exception:
        print("Consensus not available for this motif object in this Biopython version")
    # Write motif to file (JASPAR format)
    motif_file = os.path.join(OUTDIR, "demo_motif.jaspar")
    with open(motif_file, "w") as fh:
        fh.write(str(m))
    print("Motif written to:", motif_file)
except Exception as e:
    print("Motif analysis failed:", e)

print("\n=== Pipeline Finished ===")
