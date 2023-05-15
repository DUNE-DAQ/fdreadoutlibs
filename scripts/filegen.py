import numpy as np
import click
from detdataformats.wib import WIBFrame
from detdataformats.wib2 import WIB2Frame
from detdataformats.daphne import DAPHNEFrame, DAPHNEStreamFrame
from rich.console import Console
import math

console = Console()

def pattern_sin(frame, setter, max_frames):
    '''
    A sine pattern with a single frequency
    '''
    ticks = 32
    f = 62.5e6
    t = ticks / f
    freq = 552000
    x = t * np.arange(max_frames)
    y = 2000 + 1000 * np.sin(2 * math.pi * freq * x)
    for i in range(max_frames):
        fr = frame()
        for j in range(256):
            setter(fr, j, int(y[i]))
        yield fr
    return

def pattern_std(frame, setter, max_frames):
    '''
    A distribution with a predefined stdandard deviation
    '''
    ary = np.random.normal(1000, 100, (max_frames, 256))
    for i in range(max_frames):
        fr = frame()
        for j in range(256):
            setter(fr, j, int(ary[i, j]))
        yield fr
    return

def pattern_gauss(frame, setter, max_frames):
    '''
    A gaussian peak and zero everywhere else
    Every frame has the peak displaced by one sampling time compared to the previous frame
    '''
    width = 1
    amplitude = 50
    try:
        from scipy import signal
    except ModuleNotFoundError:
        print('scipy is not installed and it is needed for this pattern. Run')
        print('pip install scipy')
        exit()
    base = signal.gaussian(max_frames, width)
    base = amplitude / base.max() * base
    ary = np.repeat(base, 256).reshape((max_frames, 256))
    for i in range(ary.shape[0]):
        ary[i, :] = np.roll(ary[i], i)
        fr = frame()
        for j in range(256):
            setter(fr, j, int(ary[i, j]))
        yield fr
    return

def pattern_zero(frame, setter, max_frames):
    """
    Zero everywhere 
    """
    ary = np.zeros((max_frames, 256))
    for i in range(max_frames):
        fr = frame()
        for j in range(256):
            setter(fr, j, int(ary[i, j]))
        yield fr
    return


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--frontend', default='wib', type=click.Choice(['wib', 'wib2', 'daphne', 'daphne_stream']), help='Which frontend to use for the frames file')
@click.option('--num-frames', default=100, help='')
@click.option('--pattern', default='std', type=click.Choice(['std', 'sin', 'gauss', 'zero']), help='')
def cli(frontend, num_frames, pattern):

    pattern = {'std': pattern_std,
               'gauss': pattern_gauss,
               'sin': pattern_sin,
               'zero': pattern_zero
               }[pattern]
    frame = {'wib': WIBFrame,
             'wib2': WIB2Frame,
             'daphne': DAPHNEFrame,
             'daphne_stream': DAPHNEStreamFrame
             }[frontend]
    setter = {'wib': lambda fr, ch, val: fr.set_channel(ch, val),
              'wib2': lambda fr, ch, val: fr.set_adc(ch, val),
              'daphne': lambda fr, ch, val: fr.set_adc(ch, val),
              'daphne_stream': lambda fr, ch, val: fr.set_adc(ch, val),
              }[frontend]

    with open('gen_frames.bin', 'wb') as f:
        for fr in pattern(frame, setter, num_frames):
            f.write(fr.get_bytes())

if __name__ == "__main__":
    try:
        cli()
    except Exception as e:
        console.print_exception()
