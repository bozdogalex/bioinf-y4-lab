import sys
import requests
import argparse
from pathlib import Path


def fetch_tiny_fastq(output_dir: str):
    """
    Fetch a tiny FASTQ file for testing from a public source.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Small FASTQ from ENA
    url = "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/SRR000001/SRR000001.fastq.gz"
    
    print(f"Downloading from: {url}")
    
    output_file = output_path / "your_reads.fastq.gz"
    
    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        
        with open(output_file, 'wb') as f:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print(f"\rProgress: {percent:.1f}%", end='')
        
        print(f"\nDownloaded: {output_file}")
        print(f"Size: {output_file.stat().st_size / (1024*1024):.2f} MB")
        return output_file
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Fetch FASTQ file")
    parser.add_argument("--output", default="data/work/IrisDanila/lab03", 
                       help="Output directory")
    args = parser.parse_args()
    
    fetch_tiny_fastq(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())