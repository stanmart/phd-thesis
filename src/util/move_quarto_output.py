import os
import pathlib

out_dir = pathlib.Path(os.environ["QUARTO_PROJECT_OUTPUT_DIR"])
src_dir = out_dir / "src"
output_files = src_dir.glob("presentation/*.html")
output_files = [
    pathlib.Path(file) for file in os.environ["QUARTO_PROJECT_OUTPUT_FILES"].split()
]
for output_file in output_files:
    output_file.rename(out_dir / output_file.name)
(src_dir / "presentation").rmdir()
src_dir.rmdir()
