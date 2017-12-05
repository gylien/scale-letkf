srcfile_main = 'config/example/config.main.K'
srcfile_cycle = 'config/example/config.cycle'
srcfile_fcst = 'config/example/config.fcst'
srcfile_scale_pp = 'config/example/config.nml.scale_pp'
srcfile_scale_init = 'config/example/config.nml.scale_init'
srcfile_grads_boundary = 'config/example/config.nml.grads_boundary'
srcfile_scale = 'config/example/config.nml.scale'
srcfile_scale_user = 'config/example/config.nml.scale_user'
srcfile_obsope = 'config/example/config.obsope'
srcfile_letkf = 'config/example/config.letkf'

outfile_main = '../../../letkf.wiki/List-of-variables-in-config.main.md'
outfile_cycle = '../../../letkf.wiki/List-of-variables-in-config.cycle.md'
outfile_fcst = '../../../letkf.wiki/List-of-variables-in-config.fcst.md'
outfile_scale = '../../../letkf.wiki/Use-of-namelist-files-of-the-SCALE-model.md'
outfile_obsope = '../../../letkf.wiki/List-of-variables-in-config.nml.obsope.md'
outfile_letkf = '../../../letkf.wiki/List-of-variables-in-config.nml.letkf.md'


def get_multiline_comments(txt, startline=0, startpos=0):
    i = startline
    j = startpos
    res = ''
    try:
        if txt[i].find('#', j) >= 0:
            j = txt[i].find('#', j)
        while txt[i][j:].lstrip().startswith('#'):
            ires = txt[i][j:].lstrip().lstrip('#').strip()
            if ires.startswith('======'):
                return i+1, ''
            elif ires != '':
                if res == '':
                    res = ires
                else:
                    res += '<br>' + ires
            i += 1
            j = 0
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


if __name__ == '__main__':

    bash_to_md(srcfile_main, outfile_main)
    bash_to_md(srcfile_cycle, outfile_cycle)
    bash_to_md(srcfile_fcst, outfile_fcst)
