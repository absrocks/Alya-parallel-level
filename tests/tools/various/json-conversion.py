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
    # Postprocess
    f["postprocess"] = "alya2pos"
    # Exclude: special by default
    f["exclude"] = ["special"]
    # Converting commands
    postprocess_cmd = ""
    alya_cmd = ""
    if "alya2posCmd" in f:
        if f["alya2posCmd"] != "":
            postprocess_cmd=f["alya2posCmd"]
        f.pop("alya2posCmd")
    postprocess_cmd=postprocess_cmd.replace("ALYA2POS", "POST")
    if "command" in f:
        if f["command"] != "":
            alya_cmd=f["command"]
        f.pop("command")
    if "timeout" in f:
        f.pop("timeout")
    # Cleaning description
    f["description"]=f["description"].replace("\n"," ")
    # Comparisons
    comparisons_orig = f["comparisons"]
    comparisons = []
    for c in comparisons_orig:
        t = {}
        if "file" in c:
            t["file"] = c["file"]
            if "mpio" in t["file"]:
                f["postprocess"] = "mpio2txt"
        if "method" in c:
            t["method"] = c["method"]
        if "rounding" in t["method"]:
            t["method"] = "power"
        if "diff" in c:
            t["tolerance"] = c["diff"]
        if "tolerance" in c:
            t["tolerance"] = c["tolerance"]
        if "tol" in c:
            t["tolerance"] = c["tol"]
        if "rows" in c:
            if c["rows"] != "ALL":
                t["rows"] = c["rows"]
        if "cols" in c:
            if c["cols"] != "ALL":
                t["cols"] = c["cols"]
        comparisons.append(t)
    f["comparisons"] = comparisons
    # Executions
    executions_orig = f["executions"]
    executions = []
    for e in executions_orig:
        t = {}
        if "MPI-OMP" in e:
            t["mpi"] = e["MPI-OMP"][0]
            t["openmp"] = e["MPI-OMP"][1]
            if alya_cmd != "" or postprocess_cmd != "":
                t["commands"] = {}
            if alya_cmd != "":
                t["commands"]["alya"] = alya_cmd
            if postprocess_cmd != "":
                t["commands"]["postprocess"] = postprocess_cmd
        executions.append(t)
    f["executions"] = executions
    saveJson(f, json_file)
    print("OK")
except:
    print("Failed!")
