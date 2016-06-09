#!/usr/bin/env python
"""Tests that the creation of a Cygwin dll has the same contents as a Linux so."""

import sys
import os
import subprocess
import re

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../"

def test_cygwin_dll(prt):
    """Ensure Cygwin dll is correct."""
    os_slib = _mk_sharedlibs()
    os2data = _get_os2slibdata(os_slib)
    assert len(os2data['Linux']) > 1000
    assert os2data['CYGWIN'] == os2data['CYGWIN']
    prt.write("TEST PASSED: cyghtsdll.dll =~ libhts.so")

def _mk_sharedlibs():
    """Compile shared libs for Linux and Cygwin."""
    os.system("cd {ROOT}; make clean".format(ROOT=ROOT))
    cmds = [
        'cd {ROOT}; make clean-so libhts.so'.format(ROOT=ROOT),
        'cd {ROOT}; make clean-cygdll cyghtsdll.dll'.format(ROOT=ROOT)]
    for cmd in cmds:
        os.system(cmd)
    return [
        ('Linux', '{ROOT}/libhts.so'.format(ROOT=ROOT)),
        ('CYGWIN', '{ROOT}/cyghtsdll.dll'.format(ROOT=ROOT))]

def _get_os2slibdata(os_slib):
    """Use 'nm' to get the data from a shared lib."""
    os2data = {}
    for osname, slib in os_slib:
        data = set()
        cmd = "nm {SLIB}".format(SLIB=slib)
        results, err = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE).communicate()
        if err is not None:
            raise Exception("ERROR RUNNING({CMD}): {ERR}".format(CMD=cmd, ERR=err))
        for line in results.split('\n'):
            line = line.rstrip() # chomp
            if line:
                mtch = re.search(r'((\S)\s(\S+))$', line)
                if mtch:
                    data.add(mtch.group(1))
                else:
                    raise Exception("UNEXPECTED: {}".format(line))
        os2data[osname] = data
    return os2data

if __name__ == '__main__':
    test_cygwin_dll(sys.stdout)
