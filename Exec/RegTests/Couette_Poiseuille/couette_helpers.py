def read_input(fpath,defaults={},prefix=''):
    """Read parameters with specified prefix from input file"""
    prob_parm = {key:val for key,val in defaults.items()}
    with open(fpath,'r') as f:
        for line in f:
            line = line.strip().split('#')[0] # ignore comments
            if line == '':
                continue
            if line.startswith(prefix):
                #print('parse:',line)
                key, defn = line.split('=')
                key = key[len(prefix):].rstrip()
                value = []
                for val in defn.split():
                    try:
                        val = float(val)
                    except ValueError:
                        val = val.strip().strip('"').strip("'")
                    finally:
                        value.append(val)
                if len(value) == 1:
                    value = value[0]
                prob_parm[key] = value
    return prob_parm

