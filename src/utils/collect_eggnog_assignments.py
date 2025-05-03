"""Collect eggNOG annotations into two flat TSV files that describe their kegg and protein IDs"""

from glob import iglob
from pathlib import Path
from typing import cast

import polars as pl
from plumbum import cli, LocalPath


class EggNOGMapperResultWriter(cli.Application):
    """Write EggNOG-mapper results to Feature Count table"""

    _emapper_wdir: Path
    _kegg_ids_out: str = cast(
        str,
        cli.SwitchAttr(
            ["-k", "--kegg-ids-out"], str, help="KEGG ids output file name", default="eggnog-annotations-kegg-ids.tsv"
        ),
    )
    _by_protein_annotation: str = cast(
        str,
        cli.SwitchAttr(
            ["-b", "--by-protein-out"],
            str,
            help="By-protein annotation file name",
            default="eggnog-annotations.by-protein.tsv",
        ),
    )

    @cli.switch(["-e", "--emapper-results-dir"], cli.ExistingDirectory, mandatory=True)
    def set_emapper_wdir(self, path: LocalPath) -> None:
        """Path to `emapper_runner.py` top-level results directory"""
        self._emapper_wdir = Path(path)

    def main(self) -> None:
        """Execute"""
        eggnog_dir = self._emapper_wdir

        out = []
        ids = []
        categories = "eggNOG_OGs\tCOG_category\tDescription\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tPFAMs".split(
            "\t"
        )
        by_prot = []
        for file in iglob("*/*.emapper.annotations", root_dir=eggnog_dir):
            path = eggnog_dir.joinpath(file)
            fasta_id = file.split("/")[1].replace(".emapper.annotations", "")
            df = pl.read_csv(path, separator="\t", skip_rows=4, comment_prefix="##")
            annotations = {
                row[0][3:]: row[1]
                for row in df["KEGG_ko"].str.split(",").list.explode().drop_nulls().value_counts().iter_rows()
            }
            out.append(annotations)
            ids.append(fasta_id)

            for row in df.select("#query", *categories).iter_rows(named=True):
                by_prot.append({"File ID": fasta_id, **row})

        pl.DataFrame(by_prot).rename({"#query": "Protein ID"}).write_csv(self._by_protein_annotation, separator="\t")

        (pl.DataFrame(out).with_columns(pl.Series("id", ids)).fill_null(0.0).select("id", pl.exclude("id"))).write_csv(
            self._kegg_ids_out, separator="\t"
        )


def _main() -> None:
    EggNOGMapperResultWriter.run()


if __name__ == "__main__":
    _main()
