
#=============================================================

class InputTemplate(object):
    #--------------------------------------
    def __init__(self):
        self.original = []
        self.replace = []
        self.linenums = []

    #--------------------------------------
     # Reads in an input file and stores the information

    def readTemplate(self, filename):
        self.original[:] = []
        with open(filename,"r") as infile:
            for i, line in enumerate(infile):
                self.original.append(line)
                parts = line.split()
                if len(parts) > 0:
                    if parts[0] == "#Rep":
                        self.linenums.append(i)
#        print("Found %d replaceable lines" %len(self.linenums))

    #--------------------------------------
    def addReplacement(self, newstring):
        self.replace.append(newstring)

    #--------------------------------------
    def modifyReplacement(self, newstring, indx):
        self.replace[indx] = newstring

    #--------------------------------------
    def setReplacement(self, newstrings):
        assert isinstance(newstrings, list)
        self.replace.copy(newstrings)

    #--------------------------------------
    def clearReplacement(self):
#        self.replace.clear()
        self.replace[:] = []

    #--------------------------------------
    def generateNewFile(self, filename):
        outfile = open(filename, "w")
        itter = 0
        for i, line in enumerate(self.original):
#            print(i, itter)
            if i == self.linenums[itter]:
                outfile.write(self.replace[itter])
                if itter < len(self.linenums)-1:
                    itter += 1
            else:
                outfile.write(line)
        outfile.close()

    #--------------------------------------
    def generateNewList(self):
        itter = 0
        outlist = []
        for i, line in enumerate(self.original):
            if i == self.linenums[itter]:
                outlist.append(self.replace[itter])
                if itter < len(self.linenums)-1:
                    itter += 1
            else:
                outlist.append(line)
        return outlist

    #--------------------------------------
    def getstring(self):
        newList = self.generateNewList()
        for i, obj in enumerate(newList):
            newList[i] = obj.strip()
        newstr = '\n'.join(newList)
        return newstr

    #--------------------------------------
#=============================================================
if __name__ == "__main__":
    test = InputTemplate()
    test.readTemplate(filename="test1.txt")
    test.addReplacement("Look it is a test!")
    test.addReplacement("Yes it did!")
    newstr = test.getstring()
    print newstr
