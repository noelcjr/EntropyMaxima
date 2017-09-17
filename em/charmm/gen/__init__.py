from jinja2 import Environment, PackageLoader

env = Environment(
    loader=PackageLoader('em.charmm.gen','templates')
)


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
