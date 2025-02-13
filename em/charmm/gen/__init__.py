from jinja2 import Environment, PackageLoader

env = Environment(
    loader=PackageLoader('em.charmm.gen','templates')
)
# TODO gen_parameters(), strewam_nina_radii() gen_setup_one() need to have the
# paths as a variable that depends on installation directory

def gen_parameters():
    template = env.get_template('parameters.template')
    return template.render()
#          param_absolute_path="/home//nnoel/Code/EntropyMaxima/em/params/charmm27.ff/")

def stream_nina_radii():
    template = env.get_template('stream_nina_radii.template')
    return template.render()
# param_absolute_path="/home//nnoel/Code/EntropyMaxima/em/params/")

def fix_sequence(chain_info):
    template = env.get_template('fix_sequence.template')
    return template.render(
           a_sequence=chain_info)

def disulfide_patch(disu_bond):
    template = env.get_template('disulfides.template')
    return template.render(
           chain_num_chain_num=disu_bond)

def wmain(value):
    template = env.get_template('wmain.template')
    return template.render(val=value)

def write_psf(file_name):
    template = env.get_template('write_psf.template')
    return template.render(out_filename=file_name)

def write_psf_xplor(file_name):
    template = env.get_template('write_psf_xplor.template')
    return template.render(out_filename=file_name)

def read_psf(file_name):
    template = env.get_template('read_psf.template')
    return template.render(in_filename=file_name)

def write_pdb(file_name):
    template = env.get_template('write_pdb.template')
    return template.render(out_filename=file_name)

def write_ic(file_name):
    template = env.get_template('write_ic.template')
    return template.render(out_filename=file_name)

def write_crd(file_name):
    template = env.get_template("write_crd.template")
    return template.render(out_filename=file_name)

def read_crd(file_name):
    template = env.get_template("read_crd.template")
    return template.render(in_filename=file_name)

def read_crd2(file_name):
    template = env.get_template("read_crd2.template")
    return template.render(in_filename=file_name)

def build_heavy_atoms():
    template = env.get_template("build_heavy.template")
    return template.render()

def build_hydrogens():
    template = env.get_template("build_hydrogens.template")
    return template.render()

def stop():
    template = env.get_template("stop.template")
    return template.render()

def minimization1(outfile,input_list):
    template = env.get_template("minimization1.template")
    return template.render(
        OUTFILE=outfile,
        A=input_list[0][0],
        ROTA="resid "+input_list[0][1]+" : "+input_list[0][2],
        LINA="resid "+input_list[0][3]+" : "+input_list[0][4],
        CENA="resid "+input_list[0][5]+" : "+input_list[0][6],
        B = input_list[1][0],
        ROTB="resid "+input_list[1][1]+" : "+input_list[1][2],
        LINB="resid "+input_list[1][3]+" : "+input_list[1][4],
        CENB="resid "+input_list[1][5]+" : "+input_list[1][6])

def gen_setup_one():
    template = env.get_template('setup_one.inp.template')
    a_pair = {
        'id': 'A',
        'seq': "A.SEQ",
        'first': 'none',
        'last': 'none',
        'inp': "A_FIXRES.INP"
    }
    b_pair = {
        'id': 'B',
        'seq': "B.SEQ",
        'first': 'none',
        'last': 'none',
        'inp': "B_FIXRES.INP"
    }

    return template.render(
        param_absolute_path='/home/noel/Code/EntropyMaxima/em/params/charmm27.ff/',
        a_sequence=a_pair,
        b_sequence=b_pair,
        in_filename='INFILE',
        out_filename='OUTFILE')

if __name__ == '__main__':
    gen_setup_one()
