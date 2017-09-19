from jinja2 import Environment, PackageLoader

env = Environment(
    loader=PackageLoader('em.charmm.gen','templates')
)

def gen_parameters():
    template = env.get_template('parameters.template')
    return template.render(
           param_absolute_path='/code/em/params/charmm27.ff/')

def stream_nina_radii():
    template = env.get_template('stream_nina_radii.template')
    return template.render(
           param_absolute_path='/code/em/params/')

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

def out_psf(file_name):
    template = env.get_template('write_psf.template')
    return template.render(out_filename=file_name)

def out_psf_xplor(file_name):
    template = env.get_template('write_psf_xplor.template')
    return template.render(out_filename=file_name)

def out_pdb(file_name):
    template = env.get_template('write_pdb.template')
    return template.render(out_filename=file_name)

def out_ic(file_name):
    template = env.get_template('write_ic.template')
    return template.render(out_filename=file_name)

def out_crd(file_name):
    template = env.get_template("write_crd.template")
    return template.render(out_filename=file_name)

def in_crd(file_name):
    template = env.get_template("read_crd.template")
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
        param_absolute_path='/code/em/params/charmm27.ff/',
        a_sequence=a_pair,
        b_sequence=b_pair,
        in_filename='INFILE',
        out_filename='OUTFILE')

if __name__ == '__main__':
    gen_setup_one()
