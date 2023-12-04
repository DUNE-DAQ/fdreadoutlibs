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

def create_temp_directory(parent, test_dir):
    with tempfile.TemporaryDirectory(prefix=test_dir, dir=parent) as tmpdirname:
        print('created temporary directory', tmpdirname)
        return tmpdirname

def create_new_directory(dp, test_dir):
    dirname = tempfile.mkdtemp(dir=dp, prefix=test_dir)
    print('created new directory', dirname)

def create_local_directory(test_dir):
    localdir = Path.cwd().absolute() / test_dir
    localdir.mkdir()

def create_test_directory(dp,test_dir):
    dirname = dp / test_dir
    avx = dp / test_dir / "avx"
    naive = dp / test_dir / "naive"
    avx_vs_naive = dp / test_dir / "avx_vs_naive"
    ref_vs_avx = dp / test_dir / "ref_vs_avx"
    ref_vs_naive = dp / test_dir / "ref_vs_naive"
    make_directory(dirname)
    make_directory(avx)
    make_directory(naive)
    make_directory(avx_vs_naive)
    make_directory(ref_vs_avx)
    make_directory(ref_vs_naive)


def create_new_test(test_dir):
    test_folder = Path.cwd().absolute() / test_dir
    test_folder.mkdir()

def run_job(impl, dp, test_dir, num_frames_to_read, time_tick_offset, select_pattern, tpg_threshold):
    #os.environ["fr"] = "2"
    #os.environ["dir1"] = "/nfs/home/ivhristo/dune/rawhit/v411/patterns/patt_golden/patt_golden_wibethframes_32"
    test_dir = dp + "/" + test_dir
    os.environ["fr"] = num_frames_to_read
    os.environ["dir1"] = dp 
    os.environ["dir2"] = test_dir + "/" + impl
    os.environ["impl"] = impl.upper()
    os.environ["diff1"] = test_dir + "/avx"
    os.environ["diff2"] = test_dir + "/naive"
 
    commands = [
        #["bash", "-c", f"for i in {1..63}; do diff ${diff1}/TP_dump_*__${i}.txt ${diff2}/TP_dump_*__${i}.txt; done"],
        #["bash", "-c", "for i in {1..63}; do echo ${i}:${impl}; done;"],
        #["bash", "-c", "for i in {1..63}; do wibeth_tpg_algorithms_emulator -f ${dir1}/patt_golden_${i}_wibeth_output.bin -r false -a SimpleThreshold -i ${impl} -n ${fr} -t 499  --save-trigprim -s __${i} > log_${i}; done"],
        ["bash", "-c", "for i in {1.."+time_tick_offset+"}; do wibeth_tpg_algorithms_emulator -f ${dir1}/"+select_pattern+"_${i}_wibeth_output.bin -r false -a SimpleThreshold -i ${impl} -n ${fr} -t "+tpg_threshold+" --save-trigprim -s __${i} > log_${i}; done"],
        ["bash", "-c", "mv log_* TP_dump_*txt ${dir2}/."],
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


def run_commands(commands):
    for command in commands:
        try:
            subprocess.run(command, check=True, timeout=60)
            print(f"Success!")
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



def compare_results(dp, test_dir, time_tick_offset):
    test_dir = dp + "/" + test_dir
    diff1 = test_dir + "/avx/"
    diff2 = test_dir + "/naive/"
    res   = test_dir + "/avx_vs_naive/"
 
    commands = [
        #["bash", "-c", "for i in {1..63}; do diff "+''.join(list(glob(os.path.join(diff1,"TP_dump_*__1.txt"))))+" "+''.join(list(glob(os.path.join(diff2,"TP_dump_*__1.txt"))))+"; done"],
        #["bash", "-c", "for i in {1.."+time_tick_offset+"}; do fname=TP_dump_*__${i}.txt; diff "+''.join(list(glob(os.path.join(diff1,"TP_dump_*__1.txt"))))+" "+''.join(list(glob(os.path.join(diff2,"TP_dump_*__1.txt"))))+" > diff_${i}; done;"],
        ["bash", "-c", "for i in {1.."+time_tick_offset+"}; do diff "+diff1+'TP_dump_*__${i}.txt'+" "+diff2+'TP_dump_*__${i}.txt'+" > "+res+"diff_${i}; done;"],
    ]
    run_commands(commands)

def compare_to_reference(dp, test_dir, impl, select_pattern, time_tick_offset):
    test_dir = dp + "/" + test_dir
    diff1 = dp + "/"           # reference
    diff2 = test_dir + "/"+impl+"/"  # algorithms implementation
    res   = test_dir + "/ref_vs_"+impl+"/"

    commands = [
        ["bash", "-c", "for i in {1.."+time_tick_offset+"}; do diff "+diff1+select_pattern+'_${i}_wibeth_output_pedsub_hits.txt'+" "+diff2+'TP_dump_*__${i}.txt'+" > "+res+"diff_${i}; done;"],
    ]
    run_commands(commands) 


#------------------------------------------------------------------------------
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('file_path', type=click.Path(exists=True))
@click.option('-i', '--interactive', is_flag=True, help="Run interactive mode", default=False, show_default=True)
@click.option('-d', "--test-dir", type=str, help="Name of temporary test directory", default="test_0", show_default=True)
@click.option('-n', '--num-frames-to-read', type=str, help="Number of frames to read. Default: select all frames.", default=32, show_default=True)
@click.option('-o', '--time-tick-offset', type=str, help="Time tick of pattern start. Default: 1 (max: 63)", default=1, show_default=True)
@click.option('-p', '--select-pattern', type=str, help="Test Pattern (patt_golden, patt_pulse, patt_edge_square, patt_edge_left, patt_edge_right). Default: patt_golden", default="patt_golden", show_default=True)
@click.option('-t', '--tpg-threshold', type=str, help="Value of the TPG threshold. Default: 499", default="499", show_default=True)



def cli(file_path: str, interactive: bool,
        test_dir: str, 
        num_frames_to_read: str,
        time_tick_offset: str,
        select_pattern: str,
        tpg_threshold: str,
        ) -> None:
    script = Path(__file__).stem

    dp = Path(file_path)
    #print('File path %s' % (dp))

    create_test_directory(dp, test_dir)
    print("Comparing AVX and NAIVE...")
    run_job("avx", str(dp), test_dir, num_frames_to_read, time_tick_offset, select_pattern, tpg_threshold)
    run_job("naive", str(dp), test_dir, num_frames_to_read, time_tick_offset, select_pattern, tpg_threshold)
    compare_results(str(dp), test_dir, time_tick_offset)

    print("Comparing AVX and reference pattern...")
    compare_to_reference(str(dp), test_dir, "avx", select_pattern, time_tick_offset)

    print("Comparing NAIVE and reference pattern...")
    compare_to_reference(str(dp), test_dir, "naive", select_pattern, time_tick_offset)

    if interactive:
        import IPython
        IPython.embed(colors="neutral")

if __name__ == "__main__":

    cli()

