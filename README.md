# lncMapper
Tool to map lncRNA reads to a reference:

# Abschnitt 1: Vorbereitung und BLAST-Suche

import subprocess
import tempfile
import os
from datetime import datetime
from Bio import SeqIO

# Funktion, um Dateien im aktuellen Verzeichnis zu finden, die mit einem bestimmten Präfix beginnen
def find_files_with_prefix(prefix, directory='.'):
    return [filename for filename in os.listdir(directory) if filename.startswith(prefix)]

# Funktion, um IDs in einer FASTA-Datei zu kürzen
def shorten_ids(input_fasta, output_fasta, max_length=50):
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            record.id = record.id[:max_length]
            record.description = ""
            SeqIO.write(record, output_handle, "fasta")

# Automatisch die Query-Sequenzdateien, die Nanopore-Reads-Datei und die Referenzdatei finden
query_fasta_files = find_files_with_prefix("query_")
nanopore_reads_file = "allreads.fasta"
reference_files = find_files_with_prefix("Reference_")

if not query_fasta_files or not reference_files:
    print("Query-Sequenzdateien oder Referenzdatei wurde im aktuellen Verzeichnis nicht gefunden.")
    exit(1)

reference_file = reference_files[0]

# IDs in der Nanopore-Reads-Datei kürzen
shortened_nanopore_reads_file = "allreads_shortened.fasta"
shorten_ids(nanopore_reads_file, shortened_nanopore_reads_file)

# Verzeichnis für die Ausgabe erstellen
current_date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
output_folder = os.path.join(os.path.dirname(nanopore_reads_file), f"matching_reads_output_{current_date}")

# Dateinamen definieren
blast_database = os.path.join(output_folder, "nanopore_reads_db")

# Minimale Ausrichtungslänge und Sequenzidentitätsprozentsatz
min_alignment_length = 500
min_sequence_identity = 90.0  # 95%

# Ausgabeordner erstellen, falls er nicht existiert
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# BLAST-Datenbank für die Nanopore-Reads erstellen
subprocess.run(["makeblastdb", "-in", shortened_nanopore_reads_file, "-dbtype", "nucl", "-out", blast_database, "-parse_seqids"])

# Abschnitt 2: BLAST-Suche und Extraktion der passenden Reads

# Funktion, um BLAST-Suche durchzuführen und passende Reads zu extrahieren
def perform_blast_and_extract_reads(query_sequence, alignment_file, matching_reads_file):
    # Temporäre Datei erstellen, um die Query-Sequenz zu speichern
    with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_query_file:
        temp_query_filename = temp_query_file.name
        temp_query_file.write(">Query\n" + query_sequence + "\n")

    # Query-Sequenz gegen die BLAST-Datenbank ausrichten
    subprocess.run(["blastn", "-query", temp_query_filename, "-db", blast_database, "-out", alignment_file, "-outfmt", "7 qseqid sseqid qlen slen length nident pident evalue bitscore qstart qend sstart send"])

    # Temporäre Query-Datei löschen
    os.remove(temp_query_filename)

    # Ausrichtungen parsen und passende Reads in eine neue FASTA-Datei extrahieren
    matching_reads = {}

    with open(alignment_file, "r") as blast_output:
        best_alignment = {}
        for line in blast_output:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            alignment_length = int(fields[4])
            sequence_identity = float(fields[6])
            sseq_id = fields[1]
            if alignment_length >= min_alignment_length and sequence_identity >= min_sequence_identity:
                if sseq_id in best_alignment:
                    if alignment_length > best_alignment[sseq_id][0]:
                        best_alignment[sseq_id] = (alignment_length, line)
                else:
                    best_alignment[sseq_id] = (alignment_length, line)

    # Gefilterte Ausrichtungen in die Ausrichtungsdatei schreiben
    with open(alignment_file, "w") as blast_output:
        blast_output.write("#qseqid\tsseqid\tqlen\tslen\tlength\tnident\tpident\tevalue\tbitscore\tqstart\tqend\tsstart\tsend\n")
        for sseq_id, (_, line) in best_alignment.items():
            blast_output.write(line.strip() + "\n")

    # Passende Reads aus der ursprünglichen FASTA-Datei extrahieren und in eine neue FASTA-Datei speichern
    matching_reads_ids = set(best_alignment.keys())

    with open(matching_reads_file, "w") as output_fasta:
        for record in SeqIO.parse(shortened_nanopore_reads_file, "fasta"):
            if record.id in matching_reads_ids:
                SeqIO.write(record, output_fasta, "fasta")

    if os.path.getsize(matching_reads_file) == 0:
        print(f"Keine passenden Reads für {alignment_file} gefunden.")
    else:
        print(f"Passende Reads wurden in '{matching_reads_file}' gespeichert.")

# BLAST-Suche durchführen und passende Reads für alle Queries extrahieren
for query_fasta_file in query_fasta_files:
    query_record = SeqIO.read(query_fasta_file, "fasta")
    query_sequence = str(query_record.seq)
    query_code = query_fasta_file.split("_")[1].split(".")[0]
    alignment_file = os.path.join(output_folder, f"alignments_{query_code}.txt")
    matching_reads_file = os.path.join(output_folder, f"matching_reads_{query_code}.fasta")
    perform_blast_and_extract_reads(query_sequence, alignment_file, matching_reads_file)

# Abschnitt 3: Zusammenfassen der alignment.txt-Dateien und Anpassen der matching_reads-Dateien

# Debugging: Überprüfen der gefundenen alignment.txt-Dateien
alignment_files = find_files_with_prefix("alignments_", output_folder)
print(f"Gefundene alignment.txt-Dateien: {alignment_files}")

# Zusammenfassen der verschiedenen alignment txt files und Auflisten der unique sseqids mit der besten Übereinstimmung und pident
best_matches_dict = {}

# Alle alignment txt files durchgehen und die besten Übereinstimmungen finden
for alignment_file in alignment_files:
    print(f"Verarbeite Datei: {alignment_file}")  # Debugging: Datei wird verarbeitet
    with open(os.path.join(output_folder, alignment_file), "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            sseq_id = fields[1]
            pident = float(fields[6])
            query_id = alignment_file.split("_")[1].split(".")[0]  # Anpassung hier
            if sseq_id not in best_matches_dict or pident > best_matches_dict[sseq_id][1]:
                best_matches_dict[sseq_id] = (query_id, pident)

# Debugging: Überprüfen der gefundenen besten Übereinstimmungen
print(f"Beste Übereinstimmungen: {best_matches_dict}")

# Neue Datei erstellen und die besten Übereinstimmungen speichern
best_matches_file_path = os.path.join(output_folder, "best_matches.txt")
with open(best_matches_file_path, "w") as file:
    file.write("sseqid\tbest_query\tpident\n")
    for sseq_id, (query_id, pident) in best_matches_dict.items():
        file.write(f"{sseq_id}\t{query_id}\t{pident}\n")

print(f"Die besten Übereinstimmungen wurden in '{best_matches_file_path}' gespeichert.")

# Funktion, um die best_matches.txt-Datei zu laden und die besten Übereinstimmungen zu extrahieren
def load_best_matches(file):
    best_matches = {}
    with open(file, "r") as f:
        next(f)  # Überspringe die Kopfzeile
        for line in f:
            fields = line.strip().split("\t")
            sseq_id = fields[0]
            best_query = fields[1]
            best_matches[sseq_id] = best_query
    return best_matches

# Funktion, um die matching_reads-Dateien anzupassen
def adjust_matching_reads(best_matches_file, output_folder):
    best_matches = load_best_matches(best_matches_file)
    
    # Alle matching_reads-Dateien finden
    matching_reads_files = [f for f in os.listdir(output_folder) if f.startswith("matching_reads_") and f.endswith(".fasta")]
    
    for matching_reads_file in matching_reads_files:
        query_code = matching_reads_file.split("_")[2].split(".")[0]
        input_path = os.path.join(output_folder, matching_reads_file)
        output_path = os.path.join(output_folder, f"adjusted_{matching_reads_file}")
        
        # Reads filtern und in eine neue Datei schreiben
        with open(output_path, "w") as output_handle:
            for record in SeqIO.parse(input_path, "fasta"):
                if record.id in best_matches and best_matches[record.id] == query_code:
                    SeqIO.write(record, output_handle, "fasta")
        
        print(f"Die Datei '{output_path}' wurde erstellt.")

# Anpassung der matching_reads-Dateien durchführen
adjust_matching_reads(best_matches_file_path, output_folder)

# Abschnitt 4: Mapping der angepassten Reads und Aufräumen

# Funktion, um passende Reads mit minimap2 auf die Referenz abzubilden und sekundäre Ausrichtungen mit samtools zu entfernen
def map_and_filter_secondary_alignments(matching_reads_file, output_sam_file):
    mapped_reads_temp_sam = os.path.join(output_folder, "temp_mapped.sam")
    
    # Reads mit minimap2 auf die Referenz abbilden
    subprocess.run(["minimap2", "-ax", "map-ont", "--frag=no", "--secondary=no", reference_file, matching_reads_file], stdout=open(mapped_reads_temp_sam, 'w'))
    
    # Sekundäre Ausrichtungen mit samtools entfernen
    with open(output_sam_file, 'w') as filtered_sam:
        subprocess.run(["samtools", "view", "-h", "-F", "0x800", mapped_reads_temp_sam], stdout=filtered_sam)
    
    # Temporäre SAM-Datei löschen
    os.remove(mapped_reads_temp_sam)

# Reads für alle Queries abbilden und sekundäre Ausrichtungen entfernen
for query_fasta_file in query_fasta_files:
    query_code = query_fasta_file.split("_")[1].split(".")[0]
    adjusted_matching_reads_file = os.path.join(output_folder, f"adjusted_matching_reads_{query_code}.fasta")
    output_sam_file = os.path.join(output_folder, f"{query_code}_sites.sam")
    map_and_filter_secondary_alignments(adjusted_matching_reads_file, output_sam_file)

print("Mapping und Entfernen von sekundären Ausrichtungen abgeschlossen.")

# Unnötige Dateien löschen, außer den finalen SAM-Dateien und Ausrichtungsdateien
for filename in os.listdir(output_folder):
    full_path = os.path.join(output_folder, filename)
    if not (filename.endswith(".sam") or filename.endswith(".txt") or filename.endswith(".fasta")):
        os.remove(full_path)

print("Unnötige Dateien wurden gelöscht.")

# BLAST-Datenbank löschen
subprocess.run(["rm", "-rf", blast_database])

#Tool to generate the Circos Plot:

import os
from pycirclize import Circos
from collections import defaultdict
import pysam
import numpy as np

def find_samfiles():
    current_dir = os.getcwd()
    sam_files = [f for f in os.listdir(current_dir) if f.endswith(".sam")]
    if not sam_files:
        raise FileNotFoundError("Keine .sam-Dateien im aktuellen Verzeichnis gefunden.")
    return sam_files

def parse_samfile(samfile_path):
    chromosome_lengths = {}
    mapped_reads = []
    coverage = defaultdict(lambda: defaultdict(int))
    chrom_name_mapping = {}
    
    with pysam.AlignmentFile(samfile_path, "r") as samfile:
        # Chromosomeninformationen aus Header (@SQ) extrahieren
        for i, header in enumerate(samfile.header.get("SQ", [])):
            chrom = f"Chr. {chr(65 + i)}"  # Chromosomen alphabetisch benennen: Chr. A, Chr. B, usw.
            original_chrom = header["SN"]
            chrom_name_mapping[original_chrom] = chrom
            length = header["LN"]  # Chromosomenlänge
            chromosome_lengths[chrom] = length

        # Reads verarbeiten
        for read in samfile.fetch(until_eof=True):
            if not read.is_unmapped:  # Nur gemappte Reads berücksichtigen
                original_chrom = read.reference_name
                chrom = chrom_name_mapping.get(original_chrom, original_chrom)  # Chromosomennamen anpassen
                pos = read.reference_start if not read.is_reverse else read.reference_end  # Start oder Endposition je nach Strang
                mapped_reads.append((chrom, pos))
                coverage[chrom][pos // 100] += 1  # Coverage in 10,000-bp-Fenstern zählen
    
    return chromosome_lengths, mapped_reads, coverage

def read_coverage_data(file_path):
    coverage_data = defaultdict(lambda: defaultdict(int))
    with open(file_path, "r") as f:
        next(f)  # Überspringe die Header-Zeile
        for line in f:
            chrom, pos, cov = line.strip().split("\t")
            coverage_data[chrom][int(pos) // 100] = int(cov)
    return coverage_data

try:
    # Generiere Coverage-Daten
    samfile_paths = find_samfiles()
    print(f"Gefundene SAM-Dateien: {', '.join(samfile_paths)}")

    # Benutzerdefinierte Namen für SAM-Sektoren eingeben
    sam_sector_names = {}
    for samfile in samfile_paths:
        name = input(f"Geben Sie einen Namen für {samfile} ein: ")
        sam_sector_names[samfile] = name

    # Daten vorbereiten
    all_chromosome_lengths = {}
    all_mapped_reads = defaultdict(list)
    all_coverage = defaultdict(lambda: defaultdict(int))  # Kombinierte Coverage-Daten

    for samfile in samfile_paths:
        chromosome_lengths, mapped_reads, coverage = parse_samfile(samfile)
        all_chromosome_lengths.update(chromosome_lengths)
        all_mapped_reads[sam_sector_names[samfile]].extend(mapped_reads)
        for chrom in coverage:
            for pos in coverage[chrom]:
                all_coverage[chrom][pos] += coverage[chrom][pos]  # Coverage-Daten kombinieren

     # Coverage-Daten in eine TXT-Datei speichern
    coverage_output_file = "coverage_data.txt"
    with open(coverage_output_file, "w") as f:
        f.write("Chromosome\tPosition\tCoverage\n")
        for chrom, bins in all_coverage.items():
            for pos, cov in bins.items():
                if cov > 0:
                    f.write(f"{chrom}\t{pos * 100}\t{cov}\n")
    print(f"Coverage-Daten wurden in {coverage_output_file} gespeichert.")

    # Coverage-Daten aus der TXT-Datei lesen
    all_coverage = read_coverage_data(coverage_output_file)

    # Chromosomen und ihre Größen für Circos vorbereiten (Skalierung in Mb)
    combined_sectors = {chrom: round(max(float(length) / 10_000, 1), 6) for chrom, length in all_chromosome_lengths.items()}  # Mindestgröße 1 Mb
    
    # Füge Sektoren für jede SAM-Datei hinzu
    for name in sam_sector_names.values():
        combined_sectors[name] = 5  # Feste Größe für die Reads-Sektoren

    # Circos-Plot erstellen
    circos = Circos(sectors=combined_sectors, space=2)

    # Farben für Chromosomen und Reads (helle, kontrastreiche Farben)
    light_colors = ["#E6F0FA", "#FAE6F0", "#E6FAEA", "#FAF4E6", "#ECE6FA", "#FAE6E6"]
    name2color = {name: light_colors[i % len(light_colors)] for i, name in enumerate(combined_sectors)}
    read_colors = ["red", "blue", "green", "orange", "purple", "brown", "pink", "cyan"]

    # Sektoren erstellen und beschriften
    for sector in circos.sectors:
        track = sector.add_track((85, 90))  # Track-Bereich
        track.axis(fc=name2color[sector.name])  # Farbe setzen
        if sector.name not in sam_sector_names.values():
            track.text(sector.name, color="black", size=8)  # Chromosomen beschriften

    # Berechne den globalen maximalen Coverage-Wert
    global_max_coverage = max(max(bins.values(), default=0) for bins in all_coverage.values())

    # Coverage-Track hinzufügen
    coverage_data = []
    for chrom, bins in all_coverage.items():
        sector = circos.get_sector(chrom)
        if sector:
            coverage_track = sector.add_track((90, 100))  # Größerer Bereich für den Coverage-Track
            coverage_track.axis()
            # Berechne die x-Werte basierend auf der Position in Mbp
            x = [0] + [pos * 100 / 10_000 for pos in bins.keys()]  # Dummy-Position hinzufügen
            # Berechne die y-Werte basierend auf der Anzahl der Positionen pro 10.000 bp und normiere sie zur globalen maximalen Coverage
            y = [1.0] + [bins[pos] / global_max_coverage for pos in bins.keys()]  # Dummy-Wert hinzufügen
            # Debug-Ausgabe der x- und y-Werte sowie der globalen maximalen Coverage
            print(f"Chromosom: {chrom}")
            print(f"x-Werte: {x}")
            print(f"y-Werte: {y}")
            print(f"Globale maximale Coverage: {global_max_coverage}")
            # Zeichne die Balken mit größerer Breite
            coverage_track.bar(x, y, width=2.0, color="black")
            # Speichere die Coverage-Daten nur für Bins mit einem Wert
            coverage_data.extend([(chrom, pos * 100 / 10_000, bins[pos]) for pos in bins.keys() if bins[pos] > 0])  # Positionen in Mbp umrechnen
            # Setze die Labels für den Coverage-Track
            coverage_track.xticks_by_interval(50, label_size=8)  # Labels alle 50 Mb setzen

    # Links (Integrationsstellen) plotten
    mapping_data = []
    for i, (sector_name, mapped_reads) in enumerate(all_mapped_reads.items()):
        color = read_colors[i % len(read_colors)]
        for chrom, pos in mapped_reads:
            if chrom in all_chromosome_lengths:  # Sicherstellen, dass Chromosom existiert
                circos.link((sector_name, 0, 5), (chrom, max(pos / 10_000, 0.001), max((pos + 1) / 10_000, 0.002)), color=color)
                mapping_data.append(f"{sector_name}\t{chrom}\t{pos}")
            else:
                print(f"Warnung: Chromosom {chrom} wurde nicht gefunden.")

    # Mapping-Koordinaten in eine TXT-Datei speichern
    mapping_output_file = "mapping_coordinates.txt"
    with open(mapping_output_file, "w") as f:
        f.write("Read_Sector\tChromosome\tPosition\n")
        f.write("\n".join(mapping_data))
    print(f"Mapping-Koordinaten wurden in {mapping_output_file} gespeichert.")

    # Plot als Vektorgrafik speichern
    output_file = "circos_plot.svg"
    fig = circos.plotfig()
    fig.suptitle("Circos-Plot der Transgen-Integrationsorte", fontsize=16)
    fig.savefig(output_file, format="svg")
    print(f"Plot wurde erfolgreich als {output_file} gespeichert.")

except FileNotFoundError as e:
    print(e)
except Exception as e:
    print(f"Ein Fehler ist aufgetreten: {e}")
