# EPI2ME to WPS Report Workflow

This folder turns an EPI2ME/ONT plasmid sequencing export into customer-facing
WPS order folders with one PDF report per sample.

If several samples have the same `Order #` in the metadata sheet, their reports
go into the same `WPS Data_Order #...` folder. Samples with different order
numbers go into separate folders. Within each order folder, output filenames are
renumbered from `001`, `002`, `003`, etc. based on that order's samples rather
than the original barcode or worksheet row number.

The main instructions below are for Windows using Ubuntu/WSL. Run the commands
in the Ubuntu terminal, not in PowerShell, Command Prompt, or Anaconda Prompt.

## 1. Open Ubuntu and Activate Python

1. Open **Ubuntu** from the Windows Start menu.

2. Copy and paste this command block:

```bash
cd /mnt/c/Users/altab/plasmid_report
git pull
source .venv/bin/activate
```

After activation, the prompt should usually show `(.venv)` at the beginning.
While `.venv` is active, use `python`, not `python3`, to run the script:

```bash
python --version
minimap2 --version
samtools --version
```

If `source .venv/bin/activate` says the file does not exist, the environment
has not been set up in this folder yet. Set it up once with:

```bash
sudo apt update
sudo apt install python3 python3-venv python3-pip minimap2 samtools
python3 -m venv .venv
source .venv/bin/activate
python -m pip install pysam numpy matplotlib
```

## 2. Input Folder Location

The report command assumes all input run folders live under:

```text
/mnt/c/WPS data/
```

## 3. Prepare the Input Folder

Start with the EPI2ME output files for a sequencing run. Put the FASTA,
GenBank, BAM, FASTQ, metadata, and any optional MAF files in one folder under
`C:\WPS data\`. A clean layout looks like this:

```text
C:\WPS data\Run_2026_04_29\
  barcode01.final.fasta
  barcode02.final.fasta
  barcode01.annotations.gbk
  barcode02.annotations.gbk
  FBD...barcode01...bam
  barcode02\
    FBD...bam
  barcode01.final.fastq
  barcode02.final.fastq
  barcode01.assembly.maf        optional
  barcode02.assembly.maf        optional
  WPS_Working_Sheet_2026_04_29.xlsx
```

Each barcode must have:

- `barcodeXX.final.fasta`
- `barcodeXX.annotations.gbk`
- `barcodeXX.final.fastq`
- one raw reads FASTQ/FASTQ.GZ input whose filename contains the barcode, such
  as `barcode01.fastq.gz`. If the raw reads are split into numbered parts with
  the same base name, such as `barcode01_0.fastq.gz`, `barcode01_1.fastq.gz`,
  or `barcode01-1.fastq.gz`, `barcode01-2.fastq.gz`,
  the script aggregates all parts as one raw input. The raw FASTQ input must
  contain at least 100 reads and 100,000 read bases. Use `--use-bam` only when
  you intentionally want to run from a raw/unmapped `.bam` file instead.

Optional files:

- `barcodeXX.assembly.maf`

If a FASTQ is missing, the run still completes and the AB1 is generated from the
FASTA with default quality scores, but the sample is reported with a warning.
This `barcodeXX.final.fastq` file is the consensus FASTQ used for AB1 quality;
it is not preferred as the raw read input for alignment metrics when a larger
barcoded FASTQ is present. If multiple FASTQ files match a barcode, the script
uses the largest eligible raw input set and reports a warning listing the ignored
files. Numbered split files with the same base name are treated as one input set.
If no matching FASTQ looks like a real raw-read file, the run stops rather than
silently falling back to BAM.

The metadata sheet must contain the same barcode numbers as the data files. For
example, a row with `Barcode #` equal to `1` matches `barcode01`.
The barcode does not need leading zeroes in the sheet: `3`, `03`, `3.0`, and
`barcode3` all match files named `barcode03...`.

The metadata file should be named `WPS Working Sheet` or something very similar,
such as `WPS_Working_Sheet_2026_04_29.xlsx`. The input folder must contain
exactly one matching metadata `.xlsx`, `.csv`, or `.tsv` file.

## 4. Mixed-Contig Samples

Mixed or contaminated samples can produce more than one contig for the same
barcode and metadata row. The script detects those extra contigs, but report
generation is based on the primary contig only.

The clearest input style is to include an explicit contig suffix in each
contig-specific filename:

```text
C:\WPS data\Run_2026_04_29\
  barcode01.contig001.final.fasta
  barcode01.contig002.final.fasta
  barcode01.contig001.annotations.gbk
  barcode01.contig002.annotations.gbk
  FBD...barcode01...bam
  barcode01.final.fastq
  WPS_Working_Sheet_2026_04_29.xlsx
```

When explicit contig suffixes are present, `contig001` is treated as the primary
contig. If there is no `contig001`, the longest FASTA is used as the primary
contig. The FASTA, GenBank, AB1, per-base CSVs, plots, multimer values, host DNA
calculation, and PDF report are all generated from that primary contig. Secondary
contigs are not packaged into separate report files.

A single barcode-level BAM or FASTQ is reused for the primary contig.
Contig-specific BAM, FASTQ, or MAF files are also supported when their filenames
include the same primary contig suffix.

The script also has a looser fallback for files that do not include contig
suffixes. If multiple `final` FASTA files are found for the same barcode, it
prints that it thinks the sample has multiple contigs and assigns `contig001`,
`contig002`, etc. in sorted filename order. If there are also multiple
unlabelled GenBank, BAM, FASTQ, or MAF files with the same barcode, they are
paired to those inferred contigs by the same sorted order when the counts match.

Explicit contig suffixes are still safer when possible, because sorted filename
order is only a fallback.

## 5. Run Report

Example command:

```bash
python3 epi2me_to_final_package.py \
  --folder-name "Run_2026_04_29"
```

This processes every barcode found in `C:\WPS data\Run_2026_04_29\` and writes
the output to `C:\WPS data\Run_2026_04_29\output\`.

### Optional Commands

Use these options only when you need to choose specific barcodes, change the
output folder, or adjust alignment settings:

```bash
python3 epi2me_to_final_package.py \
  --folder-name "Run_2026_04_29" \
  --output-dir "/mnt/c/WPS data/Run_2026_04_29/output" \
  --barcodes 1 2 \
  --threads 4 \
  --sort-memory 1G
```

- `--folder-name` names the folder under `/mnt/c/WPS data/` containing all run input files.
- `--output-dir` is where the finished customer package will be written. If
  omitted, the output goes into `C:\WPS data\<folder-name>\output\`.
- `--barcodes 1 2` limits the run to barcode01 and barcode02. Omit this option to process every barcode found.
- `--multimer-denominator classified-reads` reports monomer/dimer/trimer/tetramer percentages only among reads that were close enough to 1x/2x/3x/4x plasmid length to classify. This is the default.
- `--multimer-denominator all-eligible-reads` includes eligible mapped reads that were not classifiable and adds a base-weighted `Unclassified` column to the multimer table.
- `--use-bam` uses BAM files for alignment, host DNA, read-length distribution, and multimer metrics even when a raw FASTQ/FASTQ.GZ file is detected.
- `--circ true` or `--circular-coverage` adds a diagnostic circular coverage run. Reads are aligned to a duplicated plasmid reference, projected back to the original coordinates, and written as `circular_projected_*` files under the work/report outputs. Omit it, or pass `--circ false`, to run the original linear-only pipeline.
- `--threads 4` makes alignment faster.
- `--sort-memory 1G` gives `samtools sort` more memory.

Multimer classification uses a 15% length tolerance around each plasmid
multiple. For example, a 5,000 bp contig treats reads near 5,000 bp as monomer
and reads near 10,000 bp as dimer. The PDF table labels stay simple
(`Monomer`, `Dimer`, etc.), but the displayed percentages are base-weighted so
they match the read-length distribution graph's `Total Bases (kb)` view. The
same primary mapped reads shown as mapped in the read-length graph are eligible
for multimer classification by full read length; there is no separate MAPQ or
aligned-fraction cutoff.

The report also checks the read-length distribution for evidence that the sample
is not a single contig. If a sizeable base-weighted read-length peak is not near
any 1x-4x multiple of the reported contig length, `Single Contig?` is reported
as `No`. Reads in those non-contig peaks are excluded from the multimer
calculation so an unrelated contig-size population does not distort the
monomer/dimer/trimer/tetramer percentages.

Host DNA % is calculated by aligning the same raw reads to the bundled
`E. Coli Genome.fna` reference. A read counts as host DNA when its E. coli
alignment covers more than 1,300 bp and more than 91% of the read length.
The reported Host DNA % is the percentage of primary reads that pass that host
classification rule. Base-weighted host percentages are also written to
`report_summary.json` for review, but they are not used as the headline PDF
value because a few very long host reads can otherwise dominate the number.
By default, the raw reads come from the largest appropriately sized barcoded
FASTQ/FASTQ.GZ input set. Numbered split FASTQs with the same base name are
aggregated before alignment and read-length metrics are calculated. Pass
`--use-bam` to force BAM input.

Coverage depth and per-base CSV values use primary, non-supplementary
alignments, matching the read set used for the reported mapped-bases summary.

It is okay to reuse the same `--output-dir`. If the same barcode is run again,
the script removes the previous files for that barcode and writes fresh ones.
Reports for other barcodes in the same order folder are left alone.
The default `output` folder is ignored during input discovery, so rerunning the
same folder will not treat generated alignment files as new input BAMs.

Most runs should not use `--keep-intermediates` or `--allow-aligned-input`.
Those are debugging/override options.

## 6. Report Calculation Details

The PDF report is generated from the primary contig FASTA, the matching GenBank
file, and the raw reads for the same barcode. By default, raw reads come from
the largest appropriately sized barcoded FASTQ/FASTQ.GZ input set. Split raw
FASTQs with the same base name and a trailing `-number` or `_number` are
aggregated as one input. A raw FASTQ input set must pass the sanity check of at
least 100 FASTQ records and at least 100,000 read bases. Passing `--use-bam`
forces the workflow to use the raw/unmapped BAM instead.

Plasmid alignment uses `minimap2 -ax map-ont` against the reported primary
contig FASTA, then `samtools` converts, sorts, and indexes the alignment BAM.
The read-length graph uses the same raw read lengths used for alignment. Reads
`<= 1,000 bp` are omitted from the displayed read-length graph. The graph is
base-weighted, so each bin shows `Total Bases (kb)` rather than read count.
Mapped and unmapped bars are separated by whether the primary read mapped to the
plasmid alignment BAM.

Per-base depth and coverage are calculated from primary, non-supplementary
alignments. Depth includes deletion-supporting reads. The low-confidence base
threshold is Q12. Coverage plots start at 0 on the y-axis. Coverage plots use
10% headroom above the highest value; the read-length graph uses 16% headroom so
the `Monomer` band label sits above the tallest bar.

When `--circ true` or `--circular-coverage` is enabled, the normal coverage files are left
unchanged and an extra diagnostic is added to the work directory. It creates a
doubled plasmid FASTA, aligns the same raw reads to that doubled reference, and
projects every covered doubled-reference position back with modulo arithmetic.
The diagnostic paths and coverage summary are also recorded in
`report_summary.json` under `circular_coverage_diagnostic`.

Multimer classification is based on full raw read length for primary mapped
reads, matching the mapped-read population shown in the read-length graph. There
is no separate MAPQ or aligned-fraction cutoff. A read is classified by:

```text
read_length / contig_length
```

The tolerance is 15% around each plasmid multiple:

- Monomer: within 15% of 1x contig length
- Dimer: within 15% of 2x contig length
- Trimer: within 15% of 3x contig length
- Tetramer: within 15% of 4x contig length

The PDF multimer table displays base-weighted percentages. In the default
`classified-reads` mode, percentages are calculated only over reads that fall
inside one of the 1x-4x windows. With `--multimer-denominator all-eligible-reads`,
eligible mapped reads outside those windows remain in the denominator and appear
as a base-weighted `Unclassified` column.

The `Single Contig?` field checks both FASTA record count and the read-length
distribution. A non-multimer peak is considered significant when it has at least
3 reads and at least 8% of displayed read bases, and when the peak is not near a
1x-4x plasmid multiple. Significant non-contig peak reads are excluded from the
multimer calculation so a contaminant-sized read population does not distort the
monomer/dimer/trimer/tetramer percentages.

Host DNA % is calculated by aligning the same raw reads to the bundled
`E. Coli Genome.fna` reference. A primary read counts as host when its E. coli
alignment covers more than 1,300 bp and more than 91% of the read length. The
headline PDF value is read-count based:

```text
host-classified primary reads / total primary reads * 100
```

Base-weighted host audit values are also written to `report_summary.json`, but
they are not used as the headline PDF value.

The generated `report_summary.json` and `package_summary.json` files keep audit
fields for these calculations, including host counts, host base totals, multimer
eligible/classified/unclassified counts, multimer thresholds, selected input
paths, and any warnings such as mixed-contig or split-FASTQ handling.

## Troubleshooting

If a barcode is skipped, check `run_summary.json`.

Common causes:

- The barcode exists in EPI2ME files but not in the WPS sheet.
- The WPS sheet has a barcode row but the matching EPI2ME files are missing.
- More than one unpaired file of the same type was found for the same barcode.
  Use `contig001`, `contig002`, etc. in filenames for mixed-contig samples.
- The BAM is already aligned instead of raw/unmapped.
- Two samples produce the same filename after cleanup.

If you intentionally need to debug intermediate alignment files, add:

```bat
--keep-intermediates
```

If you intentionally need to realign a BAM that already contains mapped reads,
add:

```bat
--allow-aligned-input
```

Use that only when you are sure the BAM is supposed to be realigned.
