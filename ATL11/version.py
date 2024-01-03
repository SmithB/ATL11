#!/usr/bin/env python3
import inspect
import pathlib
import subprocess
import warnings
from importlib import metadata

def softwareVersion():
    software_version=metadata.metadata("ATL11")["Version"]
    return software_version

def softwareDate():
    softwareDate='Nov 01 2020'
    return softwareDate

def softwareTitle():
    software_title=metadata.metadata("ATL11")["Summary"]
    return software_title

def softwareUrl():
    software_url=metadata.metadata("ATL11")["Home-Page"]
    return software_url

def identifier():
    identifier='atlas_l3b_is'
    return identifier

def series_version():
    series_version=metadata.metadata("ATL11")["Version"]
    return series_version

def get_git_revision_hash(
        refname: str = 'HEAD',
        short: bool = False
    ):
    """
    Get the ``git`` hash value for a particular reference

    Parameters
    ----------
    refname: str, default HEAD
        Symbolic reference name
    short: bool, default False
        Return the shorted hash value
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'rev-parse']
    cmd.append('--short') if short else None
    cmd.append(refname)
    # get output
    with warnings.catch_warnings():
        return str(subprocess.check_output(cmd), encoding='utf8').strip()

def get_git_current_branch():
    """
    Get the current ``git`` branch
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'branch', '--show-current']
    # get output
    with warnings.catch_warnings():
        return str(subprocess.check_output(cmd), encoding='utf8').strip()

def get_git_commit_date(
        date: str = "format:%Y/%m/%d %T",
    ):
    """
    Get the date of the latest ``git`` commit

    Parameters
    ----------
    date: str, default format:%Y/%m/%d %T"
        Date formatting string
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'log', '-1']
    cmd.append(f'--date={date}')
    cmd.append(f'--format="%ad"')
    # get output
    with warnings.catch_warnings():
        return str(subprocess.check_output(cmd), encoding='utf8').strip()
