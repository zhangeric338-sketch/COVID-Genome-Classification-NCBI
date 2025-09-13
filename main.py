import click
from src.download import download_dataset

@click.command()
@click.option("--virus", required=True, help="Virus name or NCBI taxon (e.g., 'sars-cov-2', 'monkeypox').")
@click.option("--output", default="data", help="Output directory for downloads.")
def main(virus, output):
    """
    Download viral genome datasets from NCBI.
    """
    try:
        dataset_path = download_dataset(virus, output)
        click.echo(f"✅ Download complete! Dataset extracted to: {dataset_path}")
    except Exception as e:
        click.echo(f"❌ Error: {e}")

if __name__ == "__main__":
    main()
