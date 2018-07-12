srcfile_main = 'config/example/config.main.K'
srcfile_cycle = 'config/example/config.cycle'
srcfile_fcst = 'config/example/config.fcst'
srcfile_scale_pp = 'config/example/config.nml.scale_pp'
srcfile_scale_init = 'config/example/config.nml.scale_init'
srcfile_grads_boundary = 'config/example/config.nml.grads_boundary'
srcfile_scale = 'config/example/config.nml.scale'
srcfile_scale_user = 'config/example/config.nml.scale_user'
srcfile_obsope = 'config/example/config.nml.obsope'
srcfile_letkf = 'config/example/config.nml.letkf'

srcfile_common_nml = '../common/common_nml.f90'

outfile_version = 'dad4dbf-(2017.11.27)'

outfile_main = '../../../letkf.wiki/List-of-variables-in-config.main-[{:s}].md'.format(outfile_version)
outfile_cycle = '../../../letkf.wiki/List-of-variables-in-config.cycle-[{:s}].md'.format(outfile_version)
outfile_fcst = '../../../letkf.wiki/List-of-variables-in-config.fcst-[{:s}].md'.format(outfile_version)
outfile_scale = '../../../letkf.wiki/Use-of-namelist-files-of-the-SCALE-model-[{:s}].md'.format(outfile_version)
outfile_obsope = '../../../letkf.wiki/List-of-variables-in-config.nml.obsope-[{:s}].md'.format(outfile_version)
outfile_letkf = '../../../letkf.wiki/List-of-variables-in-config.nml.letkf-[{:s}].md'.format(outfile_version)


def get_multiline_comments(txt, startline=0, startpos=0, commentor='#'):
    i = startline
    j = startpos
    res = ''
    ul = False
    try:
        if txt[i].find(commentor, j) >= 0:
            j = txt[i].find(commentor, j)
        while txt[i][j:].lstrip().startswith(commentor):
            if txt[i][j:].lstrip()[len(commentor):].startswith(commentor):
                ires = ''
            else:
                ires = txt[i][j:].lstrip()[len(commentor):].strip()
            if ires.startswith('======'):
                return i+1, ''
            elif ires != '':
                if ires.startswith('- '):
                    if not ul:
                        res += '<ul>'
                        ul = True
                    res += '<li>' + ires[2:].lstrip().replace('  ', ' &nbsp;') + '</li>'
                elif ul:
                    res += '</ul>' + ires.replace('  ', ' &nbsp;') 
                    ul = False
                elif res != '':
                    res += '<br>' + ires.replace('  ', ' &nbsp;') 
                else:
                    res = ires.replace('  ', ' &nbsp;') 
            i += 1
            j = 0
        if ul:
            res += '</ul>'
        if i == startline:
            i += 1
    except IndexError:
        pass
    return i, res


def bash_to_md(srcfile, outfile):
    with open(srcfile, 'r') as f:
        txt = f.read().splitlines()

    md = ''
    md_sect = ''
    i = 0
    while i < len(txt):
        if txt[i].startswith('#======'):
            if md_sect != '':
                md_sect = "\n| Variable | Default value | Explanation |\n| --- | --- | --- |" + md_sect
                md += "\n\n" + md_sect_title + "\n" + md_sect + "\n\n***"
                md_sect = ''
            i, comments = get_multiline_comments(txt, i+1, 0)
            if comments.strip() != '':
                md_sect_title = '### ' + comments
            else:
                md_sect_title = ''

        elif txt[i].find('=') >= 0:
            find_pos = txt[i].find('=')
            varname = txt[i][0:find_pos].lstrip()
            if varname.find(' ') == -1:
                i, comments = get_multiline_comments(txt, i, find_pos+2)
                default = ''
                find_pos = comments.find('(default:')
                if find_pos >= 0:
                    find_pos2 = comments.find(')', find_pos)
                    if find_pos2 >= 0:
                        default = comments[find_pos+9:find_pos2]
                        comments = comments[0:find_pos] + comments[find_pos2+1:]

                if comments.strip() != '':
                    md_sect += "\n| " + varname + " | " + default + " | " + comments + " |"

        else:
            i += 1

    if md_sect != '':
        md_sect = "\n| Variable | Default value | Explanation |\n| --- | --- | --- |" + md_sect
        md += "\n\n" + md_sect_title + "\n" + md_sect + "\n\n***"
        md_sect = ''

    print()
    print('================================================================================')
    print('Output to:', outfile)
    print('================================================================================')
    print(md)

    with open(outfile, 'w') as fo:
        fo.write(md)


def namelist_parse(srcfile):
    with open(srcfile, 'r') as f:
        txt = f.read().splitlines()

    parsed = []
    isec = -1
    ivar = -1
    inside_sect = False
    i = 0
    while i < len(txt):
        if txt[i].startswith('&'):
            title = txt[i][1:].strip()
            if title == '':
                continue
            isec += 1
            ivar = -1
            inside_sect = True
            parsed.append({})
            parsed[isec]['title'] = title
            parsed[isec]['vars'] = []
            i += 1

        elif inside_sect and txt[i].startswith('/'):
            inside_sect = False
            i += 1

        elif inside_sect and txt[i].strip().startswith('!--') and txt[i].strip().endswith('--'):
            varname = txt[i].strip()[3:-2]
            if varname.find(' ') == -1:
                ivar += 1
                parsed[isec]['vars'].append({})
                parsed[isec]['vars'][ivar]['name'] = varname
                parsed[isec]['vars'][ivar]['value'] = None
            i += 1

        elif inside_sect and txt[i].find('=') >= 0:
            while txt[i].rstrip().endswith('&') or (txt[i].rstrip().endswith(',') and (not txt[i+1].startswith('/')) and txt[i+1].find('=') < 0):
                if txt[i+1].lstrip().startswith('!'):
                    break
                txt[i] = txt[i].rstrip().rstrip('&') + txt[i+1].lstrip().lstrip('&')
                del txt[i+1]
            find_pos = txt[i].find('=')
            if find_pos >= 0:
                varname = txt[i][0:find_pos].strip()
                if varname.find(' ') == -1:
                    ivar += 1
                    parsed[isec]['vars'].append({})
                    parsed[isec]['vars'][ivar]['name'] = varname
                    parsed[isec]['vars'][ivar]['value'] = txt[i][find_pos+1:].strip().rstrip(',')
            i += 1

        else:
            i += 1

    return parsed


def namelist_src_parse(srcfile):
    with open(srcfile, 'r') as f:
        txt = f.read().splitlines()

    parsed = []
    isec = -1
    ivar = -1
    i = 0
    while i < len(txt):
        if txt[i].lstrip().startswith('!--- &'):
            title = txt[i].lstrip()[6:].strip()
            if title == '':
                break
            isec += 1
            ivar = -1
            parsed.append({})
            parsed[isec]['title'] = title
            i += 1
            if txt[i].lstrip().startswith('!'):
                i, comments = get_multiline_comments(txt, i, 0, commentor='!')
            else:
                comments = ''
            parsed[isec]['comment'] = comments.strip()
            parsed[isec]['vars'] = []

        elif isec >= 0 and txt[i].find('::') >= 0:
            while txt[i].rstrip().endswith('&'):
                if txt[i+1].lstrip().startswith('!'):
                    break
                txt[i] = txt[i].rstrip().rstrip('&') + txt[i+1].lstrip().lstrip('&')
                del txt[i+1]
            find_pos = txt[i].find('::')
            find_pos2 = txt[i].find('=', find_pos+2)
            if find_pos2 >= 0:
                varname_raw = txt[i][find_pos+2:find_pos2].strip()
                find_pos3 = varname_raw.find('(')
                if find_pos3 >= 0:
                    if not varname_raw.endswith(')'):
                        i += 1
                        continue
                    varname = varname_raw[0:find_pos3].strip()
                    dimension = varname_raw[find_pos3+1:-1].strip()
                else:
                    varname = varname_raw
                    dimension = None
                if varname.find(' ') >= 0:
                    i += 1
                    continue
                ivar += 1
                parsed[isec]['vars'].append({})
                parsed[isec]['vars'][ivar]['name'] = varname
                parsed[isec]['vars'][ivar]['dim'] = dimension
                find_pos5 = txt[i].find('!', find_pos2+1)
                if find_pos5 >= 0:
                    default = txt[i][find_pos2+1:find_pos5].strip()
                else:
                    default = txt[i][find_pos2+1:].strip()
                i, comments = get_multiline_comments(txt, i, find_pos2+1, commentor='!')
                find_pos = comments.find('(default:')
                if find_pos >= 0:
                    find_pos2 = comments.find(')', find_pos)
                    if find_pos2 >= 0:
                        default = comments[find_pos+9:find_pos2]
                        comments = comments[0:find_pos] + comments[find_pos2+1:]
                parsed[isec]['vars'][ivar]['default'] = default.strip()
                parsed[isec]['vars'][ivar]['comment'] = comments.strip()

        else:
            i += 1

    return parsed


def nml_to_md(nml, nml_config, outfile, prefix='', suffix=''):
    md = prefix

    nml_titles = [sect['title'] for sect in nml]
    for sect in nml_config:
        if sect['title'] == '' or sect['title'] not in nml_titles:
            print("[Warning] Section '{:s}' cannot be found in the source code.".format(sect['title']))
            continue
        nml_sect_idx = nml_titles.index(sect['title'])
        nml_sect = nml[nml_sect_idx]
        md += "\n#### &" + sect['title'] + "\n\n"
        if nml_sect['comment'] != '':
            md += nml_sect['comment'] + "\n\n"
        md += "| <sub>Variable</sub> | <sub>Default value</sub> | <sub>Explanation</sub> |\n| --- | --- | --- |\n"
        nml_var_names = [var['name'] for var in nml_sect['vars']]
        for var in sect['vars']:
            if var['name'] not in nml_var_names:
                print("[Warning] Variable '{:s}' cannot be found in the source code.".format(var['name']))
                continue
            nml_var_idx = nml_var_names.index(var['name'])
            nml_var = nml_sect['vars'][nml_var_idx]

            varname_display = var['name']
            if var['value'] is None:
                varname_display = '__***__ ' + varname_display
            if nml_var['dim'] is not None:
                varname_display += "<br>&nbsp;&nbsp;&nbsp;(" + nml_var['dim'] + ")"
            md += "| <sub>" + varname_display + "</sub> | <sub>" + nml_var['default'] + "</sub> | <sub>" + nml_var['comment'] + "</sub> |\n"
        md += "\n***\n"

    md += suffix

    print()
    print('================================================================================')
    print('Output to:', outfile)
    print('================================================================================')
    print(md)

    with open(outfile, 'w') as fo:
        fo.write(md)


def nml_scale_to_md(nml_scale_pp, nml_scale_init, nml_scale, outfile, prefix='', suffix=''):
    md = prefix

    md += """
***

### config.nml.scale_pp

```
"""
    for sect in nml_scale_pp:
        if sect['title'] != '':
            md_sect = ''
            for var in sect['vars']:
                if var['value'] is None:
                    md_sect += "!--" + var['name'] + "--\n"
            if md_sect != '':
                md += "&" + sect['title'] + "\n" + md_sect + "/\n\n"

    md = md[:-1] + "```\n"

    md += """
***

### config.nml.scale_init

```
"""
    for sect in nml_scale_init:
        if sect['title'] != '':
            md_sect = ''
            for var in sect['vars']:
                if var['value'] is None:
                    md_sect += "!--" + var['name'] + "--\n"
            if md_sect != '':
                md += "&" + sect['title'] + "\n" + md_sect + "/\n\n"

    md = md[:-1] + "```\n"

    md += """
***

### config.nml.scale

```
"""
    for sect in nml_scale:
        if sect['title'] != '':
            md_sect = ''
            for var in sect['vars']:
                if var['value'] is None:
                    md_sect += "!--" + var['name'] + "--\n"
            if md_sect != '':
                md += "&" + sect['title'] + "\n" + md_sect + "/\n\n"

    md = md[:-1] + "```\n"

    md += suffix

    print()
    print('================================================================================')
    print('Output to:', outfile)
    print('================================================================================')
    print(md)

    with open(outfile, 'w') as fo:
        fo.write(md)


if __name__ == '__main__':

    bash_to_md(srcfile_main, outfile_main)
    bash_to_md(srcfile_cycle, outfile_cycle)
    bash_to_md(srcfile_fcst, outfile_fcst)

    nml = namelist_src_parse(srcfile_common_nml)
    nml_scale_pp = namelist_parse(srcfile_scale_pp)
    nml_scale_init = namelist_parse(srcfile_scale_init)
#    nml_grads_boundary = namelist_parse(srcfile_grads_boundary)
    nml_scale = namelist_parse(srcfile_scale)
#    nml_scale_user = namelist_parse(srcfile_scale_user)
    nml_obsope = namelist_parse(srcfile_obsope)
    nml_letkf = namelist_parse(srcfile_letkf)

    prefix = """
**Note:** Variables that have three asterisks (__***__) in front of their names in the following tables are those to be automatically determined by the job scripts. They should have no assigned values and should be written in `config.nml.obsope` as:
```
!--VARIABLE--
```

***
"""
    nml_to_md(nml, nml_obsope, outfile_obsope, prefix)

    prefix = """
**Note:** Variables that have three asterisks (__***__) in front of their names in the following tables are those to be automatically determined by the job scripts. They should have no assigned values and should be written in `config.nml.letkf` as:
```
!--VARIABLE--
```

***
"""
    nml_to_md(nml, nml_letkf, outfile_letkf, prefix)

    prefix = """
The configuration files `config.nml.scale_pp`, `config.nml.scale_init`, and `config.nml.scale` are the same as the namelist files for **scale-rm_pp**, **scale-rm_init**, and **scale-rm** programs, respectively. Check the [SCALE-RM model documentation](http://r-ccs-climate.riken.jp/scale/doc/) for the explanation of the namelist variables.

However, in the **SCALE-LETKF**, some SCALE namelist variables are to be automatically determined by the job scripts. They should have no assigned values and should be written in `config.nml.scale*` as:
```
!--VARIABLE--
```
These automatically determined variables are listed as follow.
"""
    nml_scale_to_md(nml_scale_pp, nml_scale_init, nml_scale, outfile_scale, prefix)
