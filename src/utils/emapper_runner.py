#!/usr/bin/env python
"""Run EggNOG mapper annotation protocol"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from glob import iglob
from pathlib import Path
from typing import cast

from plumbum import LocalPath, cli, local
from tqdm import tqdm


def emapper(fasta_file: Path, wdir: Path, emapper_data_dir: Path) -> None:
    """
    Run `emapper.py` on FASTA file

    :param fasta_file: Path to FASTA file
    :param wdir: Working directory
    :param emapper_data_dir: Path to EggNOG data
    :return: None
    """
    wdir = wdir.joinpath(fasta_file.stem)
    out_file = wdir.joinpath(f"{fasta_file.stem}.emapper.annotations")
    if not wdir.exists():
        wdir.mkdir(parents=True)
    # If already done, skip
    done_file = wdir.joinpath("done")
    if done_file.exists() or out_file.exists():
        return
    local["emapper.py"][
        "--cpu",
        "8",
        "--override",
        "-i",
        fasta_file,
        "--itype",
        "proteins",
        "-m",
        "diamond",
        "--evalue",
        "0.001",
        "--query_cover",
        "80",
        "--subject_cover",
        "80",
        "--sensmode",
        "mid-sensitive",
        "--data_dir",
        emapper_data_dir,
        "--output_dir",
        wdir,
        "--output",
        fasta_file.stem,
    ]()
    with open(done_file, "w") as done_ptr:
        done_ptr.write("done\n")


class EMapperRunner(cli.Application):
    """Annotate genomes with eggNOG mapper"""

    input_dir: Path
    wdir: Path
    emapper_data_dir: Path
    file_glob: str = cast(str, cli.SwitchAttr(["-f", "--file-glob"], str, default="*.faa", help="Input file glob"))
    n_workers: int = cast(
        int,
        cli.SwitchAttr(
            ["-n", "--n-workers"],
            cast(type, cli.Range(1, 256)),
            default=25,
            help="Number of workers (10 threads per worker)",
        ),
    )

    @cli.switch(["-i", "--input-dir"], cli.ExistingDirectory, mandatory=True)
    def set_input_dir(self, path: LocalPath):
        """Path to input directory"""
        self.input_dir = Path(path)

    @cli.switch(["-w", "--wdir"], LocalPath, mandatory=True)
    def set_wdir(self, path: LocalPath):
        """Working directory"""
        self.wdir = Path(path)

    @cli.switch(["-e", "--emapper-data"], cli.ExistingDirectory, mandatory=True)
    def set_emapper_data_dir(self, path: LocalPath):
        """Path to eggnog-mapper-data directory"""
        self.emapper_data_dir = Path(path)

    def main(self) -> None:
        """Execute"""
        if not self.wdir.exists():
            self.wdir.mkdir(parents=True)
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = []
            for genome in iglob(self.file_glob, root_dir=self.input_dir):
                fasta_file = self.input_dir.joinpath(genome)
                futures.append(executor.submit(emapper, fasta_file, self.wdir, self.emapper_data_dir))
            for result in tqdm(as_completed(futures), total=len(futures)):
                if (exc := result.exception()) is not None:
                    print(exc)


if __name__ == "__main__":
    EMapperRunner.run()
