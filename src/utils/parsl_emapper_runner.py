"""Run eggNOG-mapper on set of genomes"""

import shutil
from glob import iglob
from pathlib import Path
from typing import Iterator, cast

import parsl
import polars as pl
from parsl import python_app
from parsl.channels import LocalChannel
from parsl.launchers import SingleNodeLauncher, SrunLauncher
from parsl.providers import LocalProvider

### Data and code dependencies

INPUT_DIR: Path = Path("/pscratch/sd/c/cjneely/bacdive-kg-microbe/archaea_bacteria")
assert INPUT_DIR.is_dir()

EXISTING_ANNOTATIONS_FILE: Path = Path(
    "/global/cfs/cdirs/kbase/ke_prototype/traits/bacdive-annotations-8501/all-annotation-results.merged.tsv"
)
assert EXISTING_ANNOTATIONS_FILE.is_file()

EGGNOG_MAPPER_DATA: Path = Path("/pscratch/sd/c/cjneely/data/eggnog")
assert EGGNOG_MAPPER_DATA.is_dir()

EMAPPER_RUNNER = Path(__file__).parent.joinpath("emapper_runner.py")
assert EMAPPER_RUNNER.is_file()
EMAPPER_COLLECT = Path(__file__).parent.joinpath("collect_eggnog_assignments.py")
assert EMAPPER_COLLECT.is_file()

### Config

N_BATCHES: int = 5  # One batch per node
WDIR: Path = Path("/pscratch/sd/c/cjneely/bacdive-kg-microbe/wdir")
WDIR.mkdir(parents=True, exist_ok=True)


class _Batcher:
    def __init__(self, wdir: Path, n_batches: int):
        self._wdir_root: Path = wdir
        self._n_batches: int = n_batches
        self._it_pos: int = -1
        self._create_batch_dirs()

    def _batch_dir(self, i: int) -> Path:
        return self._wdir_root.joinpath(f"batch-dir-{i}")

    def _wdir(self, i: int) -> Path:
        return self._wdir_root.joinpath(f"emapper-wdir-{i}")

    def _create_batch_dirs(self) -> None:
        for i in range(self._n_batches):
            self._batch_dir(i).mkdir(parents=True, exist_ok=True)

    def __iter__(self) -> Iterator:
        self._it_pos = -1
        return self

    def __next__(self) -> Path:
        self._it_pos += 1
        self._it_pos %= self._n_batches
        return self._batch_dir(self._it_pos)

    def directories(self) -> Iterator[tuple[Path, Path]]:
        for i in range(self._n_batches):
            yield self._batch_dir(i), self._wdir(i)  ### Begin


parsl.load(
    parsl.Config(
        executors=[
            parsl.HighThroughputExecutor(
                label="BacDive-eggNOG_headless",
                max_workers=1,
                provider=LocalProvider(
                    channel=LocalChannel(script_dir="."),
                    nodes_per_block=N_BATCHES,
                    launcher=cast(SingleNodeLauncher, SrunLauncher(overrides="-c 128")),
                    worker_init="conda activate kg_microbe-eggnog-mapper-env",
                    cmd_timeout=120,
                    init_blocks=1,
                    max_blocks=1,
                ),
            )
        ],
        strategy=None,
    )
)


def _split_input_directory(
    wdir: Path, n_batches: int, input_dir: Path, existing_annotations_file: Path
) -> list[tuple[Path, Path]]:
    idx = "File ID"
    all_genome_ids = {
        file_id.removesuffix(".faa"): file_id.split("_")[1] for file_id in iglob("*.faa", root_dir=input_dir)
    }
    print("Existing annotations", pl.read_csv(existing_annotations_file, separator="\t"))
    completed_ids = pl.scan_csv(existing_annotations_file, separator="\t").select(idx).unique().collect()[idx]
    needs_annotation = pl.DataFrame({"basename": all_genome_ids.keys(), "gca/f": all_genome_ids.values()}).filter(
        ~pl.col("gca/f").str.split("_").list.get(1).is_in(completed_ids)
    )
    print("Needs annotation", needs_annotation)
    batcher = _Batcher(wdir, n_batches)
    for file_id, batch_wdir in zip(needs_annotation["basename"], batcher, strict=False):
        shutil.copy(input_dir.joinpath(f"{file_id}.faa"), batch_wdir)
    return list(batcher.directories())


@python_app
def _run_emapper_on_dir(
    input_dir: Path, wdir: Path, eggnog_data_path: Path, runner_path: Path, collect_path: Path
) -> Path:
    from glob import glob
    from pathlib import Path

    from plumbum import local

    local[runner_path]["-i", input_dir, "-f", "*.faa", "-w", wdir, "-n", "16", "-e", eggnog_data_path]()
    local[collect_path][
        "-e", wdir, "-k" f"{wdir}-eggnog-annotations-kegg-ids.tsv", "-b", f"{wdir}-eggnog-annotations.by-protein.tsv"
    ]()
    output_file: list[str] = glob("*.emapper.annotations", root_dir=wdir)
    assert len(output_file) == 1
    return Path(wdir.joinpath(output_file[0]))


if __name__ == "__main__":
    futures = []
    for batch_dir, wdir in _split_input_directory(WDIR, N_BATCHES, INPUT_DIR, EXISTING_ANNOTATIONS_FILE):
        futures.append(_run_emapper_on_dir(batch_dir, wdir, EGGNOG_MAPPER_DATA, EMAPPER_RUNNER, EMAPPER_COLLECT))
    for future in futures:
        future.result()
