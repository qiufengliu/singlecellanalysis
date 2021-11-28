#!/usr/bin/env python3
import click
from pathlib import Path

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("--fq1s", type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), multiple=True,
              required=True)
@click.option("--fq2s", type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), multiple=True,
              required=True)
@click.option("--sample_name", type=click.STRING, required=True)
def cli(fq1s: list, fq2s: list, sample_name: str):
    """
    change raw fastq file name to tenx file name
    """
    fq1s = [Path(i) for i in fq1s]
    fq2s = [Path(i) for i in fq2s]
    Path("rawdata").mkdir(parents=True, exist_ok=True)
    for count, (fq1, fq2) in enumerate(zip(sorted(fq1s), sorted(fq2s))):
        fq1.link_to(f"rawdata/{sample_name}_S1_L{count + 1:03d}_R1_001.fastq.gz")
        fq2.link_to(f"rawdata/{sample_name}_S1_L{count + 1:03d}_R2_001.fastq.gz")


if __name__ == "__main__":
    cli()
