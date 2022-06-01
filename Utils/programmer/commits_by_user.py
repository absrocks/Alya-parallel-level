# -*- coding: utf-8 -*-

# ==== IMPORTS SECTION ========================================================
import re
import os
import copy

# ==== CONSTANTS DEFINITIONS ==================================================

# ==== CLASS DEFINITION =======================================================


class CoomitsByUser():

    # _________________________________________________________________________
    def __init__(self):

        self.user_regex = "bsc\d{5}"
        self.users = {}

        print("Getting svn log information, this will take few minuts")
        self.svnlog = os.popen('svn log -v ../../').read()

        # Read data
        # with open("svnlog") as fp:
        #    self.svnlog = fp.read()

        # svn log sepator between commits
        self.separator = "-"*72
        self.modules = {
                "alefor": 0,
                "chemic": 0,
                "exmedi": 0,
                "fakemo": 0,
                "gotita": 0,
                "helmoz": 0,
                "immbou": 0,
                "insitu": 0,
                "levels": 0,
                "magnet": 0,
                "nastal": 0,
                "nastin": 0,
                "neutro": 0,
                "partis": 0,
                "porous": 0,
                "quanty": 0,
                "radiat": 0,
                "solidz": 0,
                "temper": 0,
                "turbul": 0,
                "wavequ": 0,
                "kernel": 0
            }

    # _________________________________________________________________________
    def split(self):
        self.svnlog = self.svnlog.split(self.separator)

    # _________________________________________________________________________
    def get_user(self, commit):
        try:
            return re.search(self.user_regex, commit).group(0)
        except Exception:
            return None

    # _________________________________________________________________________
    def get_modules(self, commit):
        for module in self.modules:
            if module in commit:
                yield module

    # _________________________________________________________________________
    def extract(self):
        self.split()
        for commit in self.svnlog:

            # Get user
            user = self.get_user(commit)
            if user is None:
                continue

            # Create user structure if user do not exist
            if user not in self.users:
                self.users[user] = {}
                self.users[user]["modules"] = copy.deepcopy(self.modules)

            # get user modules
            user_modules = self.users[user]["modules"]
            commited_modules = self.get_modules(commit)

            for module in commited_modules:
                user_modules[module] += 1

    # _________________________________________________________________________
    def report(self):
        for user in self.users:
            print("User: " + user)
            modules = self.users[user]["modules"]
            for module in modules:
                if modules[module] > 0:
                    print("    - " + module + ": " + str(modules[module]))
            print()

    # _________________________________________________________________________
    def generate_csv(self):
        name = "commits_by_user.csv"
        headers = ["User"]

        headers.extend([module for module in sorted(self.modules)])
        fp = open(name, "w")
        fp.write(",".join(headers)+"\n")

        for user in self.users:
            modules = self.users[user]["modules"]
            fp.write(user+",")
            fp.write(
                ",".join([str(modules[module]) for module in sorted(modules)]))
            fp.write("\n")


if __name__ == "__main__":
    commits = CoomitsByUser()
    commits.extract()
    commits.report()
    commits.generate_csv()
    print("csv commits_by_user.csv generated!")
