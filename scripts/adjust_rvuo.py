import gzip
import sys
import re

def get_contig_offsets_and_lengths(contig_file):
    offsets = {}
    chrom_lengths = {}
    try:
        with open(contig_file, 'r') as f:
            for line in f:
                if line.strip():
                    contig_name, start, end = line.strip().split('\t')
                    start = int(start)
                    end = int(end)
                    # Extract chromosome (e.g., '9' from SWPC#0#9)
                    chrom_match = re.match(r'RVUO#0#(\d+)', contig_name)
                    if not chrom_match:
                        print(f"Warning: Invalid contig name format: {contig_name}")
                        continue
                    chrom = chrom_match.group(1)
                    # Update max end position for this chromosome
                    chrom_lengths[chrom] = max(chrom_lengths.get(chrom, 0), end)
                    # Assign piece number for this segment
                    segment_count = len([k for k in offsets if k.startswith(contig_name + '#piece')])
                    piece_name = f"{contig_name}#piece{segment_count}"
                    offsets[piece_name] = start - 1  # Store offset (0-based)
        return offsets, chrom_lengths
    except Exception as e:
        print(f"Error reading contig file: {e}")
        sys.exit(1)

def extract_chrom_from_contig(contig):
    # Extract chromosome from contig name like SWPC#0#9#piece0
    try:
        parts = contig.split('#')
        if len(parts) >= 3:
            return parts[2]  # e.g., '9' from SWPC#0#9#piece0
        else:
            raise ValueError(f"Invalid contig format: {contig}")
    except Exception as e:
        print(f"Error extracting chromosome from contig {contig}: {e}")
        return None

def adjust_vcf(input_vcf, output_vcf, contig_file, compress_output=False):
    offsets, chrom_lengths = get_contig_offsets_and_lengths(contig_file)
    try:
        # Open input file
        if input_vcf.endswith('.gz'):
            infile = gzip.open(input_vcf, 'rt')
        else:
            infile = open(input_vcf, 'r')

        # Open output file
        if compress_output or output_vcf.endswith('.gz'):
            outfile = gzip.open(output_vcf, 'wt')
        else:
            outfile = open(output_vcf, 'w')

        with infile, outfile:
            for line in infile:
                if line.startswith('##contig='):
                    # Parse contig header
                    contig_match = re.match(r'##contig=<ID=([^,]+),length=\d+>', line.strip())
                    if contig_match:
                        contig_id = contig_match.group(1)
                        chrom = extract_chrom_from_contig(contig_id)
                        if chrom and chrom in chrom_lengths:
                            # Write new contig header with chromosome and max length
                            outfile.write(f"##contig=<ID={chrom},length={chrom_lengths[chrom]}>\n")
                        continue  # Skip original contig line
                    else:
                        print(f"Warning: Invalid contig header format: {line.strip()}")
                        continue
                elif line.startswith('#'):
                    # Preserve other header lines
                    outfile.write(line)
                    continue
                
                fields = line.strip().split('\t')
                contig = fields[0]
                try:
                    pos = int(fields[1])
                except ValueError:
                    print(f"Warning: Invalid position in line: {line.strip()}")
                    continue

                # Extract chromosome from contig name
                original_chrom = extract_chrom_from_contig(contig)
                if original_chrom is None:
                    print(f"Warning: Skipping line with invalid contig {contig}: {line.strip()}")
                    continue

                # Find matching contig piece
                if contig in offsets:
                    offset = offsets[contig]
                    new_pos = pos + offset
                    fields[0] = original_chrom
                    fields[1] = str(new_pos)
                    outfile.write('\t'.join(fields) + '\n')
                else:
                    print(f"Warning: Unknown contig {contig}, skipping line: {line.strip()}")

    except Exception as e:
        print(f"Error processing VCF: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python adjust_vcf.py input_vcf output_vcf contig_file compress_output")
        print("Example: python adjust_vcf.py input.vcf.gz output.vcf.gz contigs.txt true")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    contig_file = sys.argv[3]
    compress_output = sys.argv[4].lower() == "true"

    adjust_vcf(input_vcf, output_vcf, contig_file, compress_output)