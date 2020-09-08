import subprocess

output = subprocess.getoutput("qstat -u mcanyelles")

if output == '':
    string = 'No current jobs'

else :
    output = output.split('\n')[5:]
    
    string = ''

    for line in output:
        string = string + '\n' + \
                    line[:6] + ' | ' + \
                    line[30:38] + ' | ' + \
                    line[39:56] + ' | ' + \
                    line[86:87] + ' | ' + \
                    line[88:]
    rs = int(string.count('| R |'))
    string = '%s calcs are running:\n' % rs + string

print(string)
