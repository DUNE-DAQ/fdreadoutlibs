import click
from pathlib import Path
import subprocess
import tempfile
from glob import glob

import os

def make_directory(dirname):
    try:
        dirname.mkdir()
    except FileExistsError as error:
        print(error)

def create_temp_directory(parent, name):
    with tempfile.TemporaryDirectory(prefix=name, dir=parent) as tmpdirname:
        print('created temporary directory', tmpdirname)
        return tmpdirname

def create_new_directory(dp, name):
    dirname = tempfile.mkdtemp(dir=dp, prefix=name)
    print('created new directory', dirname)

def create_local_directory(name):
    localdir = Path.cwd().absolute() / name
    localdir.mkdir()

def create_test_directory(dp,name):
    dirname = dp / name
    avx = dp / name / "avx"
    naive = dp / name / "naive"
    make_directory(dirname)
    make_directory(avx)
    make_directory(naive)


def create_new_test(name):
    test_folder = Path.cwd().absolute() / name
    test_folder.mkdir()

def run_job(impl, name):
    os.environ["fr"] = "2"
    os.environ["dir1"] = "/nfs/home/ivhristo/dune/rawhit/v411/patterns/patt_golden/patt_golden_wibethframes_32"
    os.environ["dir2"] = name + "/" + impl
    os.environ["impl"] = impl.upper()
    os.environ["diff1"] = name + "/avx"
    os.environ["diff2"] = name + "/naive"
 
    commands = [
        #["bash", "-c", f"for i in {1..63}; do diff ${diff1}/TP_dump_*__${i}.txt ${diff2}/TP_dump_*__${i}.txt; done"],
        #["bash", "-c", "for i in {1..63}; do echo ${i}:${impl}; done;"],
        ["bash", "-c", "for i in {1..63}; do wibeth_tpg_algorithms_emulator -f ${dir1}/patt_golden_${i}_wibeth_output.bin -r false -a SimpleThreshold -i ${impl} -n ${fr} -t 499  --save-trigprim -s __${i} > dbg_${i}; done"],
        ["bash", "-c", "mv dbg_? dbg_?? TP_dump_*txt ${dir2}/."],
    ]
    for command in commands:
        try:
            subprocess.run(command, check=True, timeout=60)
        except FileNotFoundError as exc:
            print(
                f"Command {command} failed because the process "
                f"could not be found.\n{exc}"
            )
        except subprocess.CalledProcessError as exc:
            print(
                f"Command {command} failed because the process "
                f"did not return a successful return code.\n{exc}"
            )
        except subprocess.TimeoutExpired as exc:
            print(f"Command {command} timed out.\n {exc}")

def compare_results(name):
    diff1 = name + "/avx/"
    diff2 = name + "/naive/"
 
    commands = [
        ["bash", "-c", "for i in {1..63}; do diff "+''.join(list(glob(os.path.join(diff1,"TP_dump_*__1.txt"))))+" "+''.join(list(glob(os.path.join(diff2,"TP_dump_*__1.txt"))))+"; done"],
    ]
    for command in commands:
        try:
            subprocess.run(command, check=True, timeout=60)
        except FileNotFoundError as exc:
            print(
                f"Command {command} failed because the process "
                f"could not be found.\n{exc}"
            )
        except subprocess.CalledProcessError as exc:
            print(
                f"Command {command} failed because the process "
                f"did not return a successful return code.\n{exc}"
            )
        except subprocess.TimeoutExpired as exc:
            print(f"Command {command} timed out.\n {exc}")



#------------------------------------------------------------------------------
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('file_path', type=click.Path(exists=True))
@click.option('-i', '--interactive', is_flag=True, help="Run interactive mode", default=False, show_default=True)
@click.option('-n', "--name", type=str, help="Name of temporary test directory", default="test_0", show_default=True)

def cli(file_path: str, interactive: bool,
        name: str, 
        ) -> None:
    script = Path(__file__).stem

    dp = Path(file_path)
    print('File path %s' % (dp))

    create_test_directory(dp, name)
    run_job("avx", name)
    run_job("naive", name)
    compare_results(name)

    if interactive:
        import IPython
        IPython.embed(colors="neutral")

if __name__ == "__main__":

    cli()

