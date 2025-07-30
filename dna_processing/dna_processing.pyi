from typing import Set, List, Tuple

def encode_fasta(
    chromosomes: Set[str],
    fasta_file: str,
    output_directory: str
) -> None: ...

def mutate(
    chromosome: str,
    mutations: List[Tuple[int, str, str]],
    chromosome_data_directory: str,
    output_filepath: str
) -> bool: ...
