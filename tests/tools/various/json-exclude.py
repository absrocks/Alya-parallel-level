#Author: Damien Dosimont

import sys
import os
import json

def openJson(file):
    try:
        with open(file) as json_data:
            f = json.load(json_data)
    except IOError:
        print("Could not open the JSON config file: "+file)
        raise
    except ValueError:
        print("File format Error in the JSON config file: "+file)
        raise
    except Exception:
        print("Error while opening the JSON file: "+file)
        raise
    else:
        return f

def saveJson(d, name):
    """
        Save a json file
    """
    if(type(d) != type({})):
        logger.error("Date type is not a dict")
        raise TypeError("First parameter data type expected: dict " +
                        "but find: "+type(d))
    try:
        with open(name, 'w') as fp:
            json.dump(d, fp, indent=4, sort_keys=True)
    except IOError:
        logger.error("Could not open the JSON config file: "+name)
        raise
    except Exception:
        logger.error("Error while opening the JSON file: "+name)
        raise
    else:
        return 1

json_file = sys.argv[1]
try:
    f = openJson(json_file)
    f["exclude"] = ["special"]
    saveJson(f, json_file)
    print("OK")
except:
    print("Failed!")
