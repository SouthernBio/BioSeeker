from os import system, name


def delete_csv_files():
    """Deletes unnecessary CSV files from the current working directory"""
    if name == 'posix':
        system('rm history_*.csv')
        system('clear')
    else:
        system('del history_*.csv')
        system('cls')

    system('echo ALL FILES ASSEMBLED.')

    return