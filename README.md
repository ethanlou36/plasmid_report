# EPI2ME to WPS Report Workflow

This folder turns an EPI2ME/ONT plasmid sequencing export into customer-facing
WPS order folders with one PDF report per sample.

If several samples have the same `Order #` in the metadata sheet, their reports
go into the same `WPS Data_Order #...` folder. Samples with different order
numbers go into separate folders.

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
- `--multimer-denominator all-eligible-reads` includes eligible mapped reads that were not classifiable and adds an `Unclassified` column to the multimer table.
- `--use-bam` uses BAM files for alignment, host DNA, read-length distribution, and multimer metrics even when a raw FASTQ/FASTQ.GZ file is detected.
- `--threads 4` makes alignment faster.
- `--sort-memory 1G` gives `samtools sort` more memory.

Multimer classification uses a 15% length tolerance around each plasmid
multiple. For example, a 5,000 bp contig treats reads near 5,000 bp as monomer
and reads near 10,000 bp as dimer. The PDF table labels stay simple
(`Monomer`, `Dimer`, etc.), but the displayed percentages are base-weighted so
they match the read-length distribution graph's `Total Bases (kb)` view.

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

It is okay to reuse the same `--output-dir`. If the same barcode is run again,
the script removes the previous files for that barcode and writes fresh ones.
Reports for other barcodes in the same order folder are left alone.
The default `output` folder is ignored during input discovery, so rerunning the
same folder will not treat generated alignment files as new input BAMs.

Most runs should not use `--keep-intermediates` or `--allow-aligned-input`.
Those are debugging/override options.

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
