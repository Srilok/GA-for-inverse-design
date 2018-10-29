#import subprocess as subproc
import subprocess32 as subproc
import os
import time
import sys 
#import re
#from mpi4py import MPI


#==================================================================================
from contextlib import contextmanager

@contextmanager
def stdout_redirector(stream):
    old_stdout = sys.stdout
    sys.stdout = stream
    try:
        yield
    finally:
        sys.stdout = old_stdout

#==================================================================================
class Simulation(object):
    def __init__(self, name="Program", sdtin=subproc.PIPE, stdout = subproc.PIPE):
        self.template = []
        self._progArguments = []
        self.progDirectory = ""
        self.name = name
        self.stdin = sdtin
        self.stdout = stdout

    def addProgramArg(self, command):
        self._progArguments.append(command)

    def setProgArg(self, commandList):
        self._progArguments = commandList

    def clearProgArg(self):
#        self._progArguments.clear()
        self._progArguments[:] = []

    def getProgArg(self):
        return self._progArguments

    def setWorkDir(self, directory):
        self.progDirectory = directory

    def getWorkDir(self):
        return self.progDirectory

    def loadForceTemplate(self, fileName):
        with open(fileName) as infile:
            self.template = infile.readlines()

    def changeStdIn(self, file):
        self.stdin = file

    def runProgram(self, echoOutput=False, redirect_file = sys.stdout):

        with stdout_redirector(open(redirect_file,"a")): 
            print("Running Command:", self._progArguments)
            try:
                prog = subproc.Popen(self._progArguments, stdin=self.stdin, stdout=subproc.PIPE, stderr=subproc.PIPE,
                                     shell=False, universal_newlines=True,
                                     cwd=self.progDirectory, preexec_fn=os.setsid)
                if echoOutput:
                    for line in prog.stdout:
                        print(line)
                rtnCode = prog.wait()
                if rtnCode:
                   raise subproc.CalledProcessError(rtnCode, self._progArguments)
            except:
                print(self.name + " did not exit properly!")
                raise
                prog.terminate()
                return -1
            else:
                print(self.name + " exited with return code 0.")
                prog.terminate()
                return 0

    def runProgram_MPI(self, echoOutput=False, nProc=1, timeout = -1 ):
#        mpiArgs = ["/opt/openmpi/bin/mpirun"]
        mpiArgs = ["mpirun"]
        mpiArgs.append("-np")
        mpiArgs.append( str(nProc) )
        try: 
            nodeFile = os.environ['PBS_NODEFILE'] 
            mpiArgs.append("-machinefile")
            mpiArgs.append(str(nodeFile))
        except:
            nodeFile = None
        mpiArgs = mpiArgs + self._progArguments
        print("Running Command:", mpiArgs)
        try:
            prog = subproc.Popen(mpiArgs, stdin=self.stdin, stdout=subproc.PIPE, stderr=subproc.PIPE,
                                 shell=False, universal_newlines=True,  cwd=self.progDirectory, preexec_fn=os.setsid)
            if echoOutput:
                for line in prog.stdout:
                    print(line)
            if timeout > 0:
                print "Time Out Value:", timeout
                try:
#                    stdout, stderr = prog.check_output(timeout=timeout)
                    rtnCode = prog.communicate(timeout=timeout)
                except subproc.TimeoutExpired:
                    print self.name + " timed out!"
                    prog.kill()
            else:
                rtnCode = prog.communicate()

            if echoOutput:
               for line in prog.stdout:
                   print line

            if rtnCode > 0:
                raise subproc.CalledProcessError(rtnCode, self._progArguments)
            prog.terminate()
        except:
            print(self.name + " did not exit properly!")
            prog.terminate()
            return -1
        else:
            print(self.name + " exited with return code 0.")
            prog.terminate()
            return 0

    def runProgram_OMP(self, echoOutput=False, nProc=1):
        print("Running Command:", self._progArguments)
        ompEnv = dict(os.environ)
        ompEnv["OMP_NUM_THREADS"] = str(nProc)
        try:
            prog = subproc.Popen(self._progArguments, stdin=self.stdin, stdout=subproc.PIPE, shell=False, universal_newlines=True,
                                 cwd=self.progDirectory, env=ompEnv)
            if echoOutput:
                for line in prog.stdout:
                    print(line)
            rtnCode = prog.wait()
            if rtnCode > 0:
                raise subproc.CalledProcessError(rtnCode, self._progArguments)
        except:
            print(self.name + " did not exit properly!")
            prog.terminate()
            return -1
        else:
            print(self.name + " exited with return code 0.")
            prog.terminate()
            return 0
#==================================================================================
class LAMMPS_Sim(Simulation):
    def __init__(self):
        Simulation.__init__(self, name="LAMMPS")
#        self.addProgramArg("lmp_mpi")
        self.addProgramArg("/home/share/cnm50256/bin/lammps-31Mar17-mod/bin/lmp_mpi_3")
        self.__inputScript = subproc.PIPE

    def setInput(self, inputName):
        self.__inputScript = open(inputName, "r")
        self.changeStdIn(self.__inputScript)

    def closeInput(self):
        self.__inputScript.close()

#==================================================================================
class Nucleation_Sim(Simulation):
    def __init__(self):
        Simulation.__init__(self, name="Nucleation")
        self.setWorkDir(os.getcwd())
        self.exeDir = "/home/share/cnm50256/WTworkflow/scripts/exec/"
#        self.exeDir = "/Users/tloeffler/SimulationCodes/GitNucleation/"
        self.addProgramArg(self.exeDir+"generalNucleation")
        self.inputScript = "dummy.dat"
#        self.addProgramArg("TersofInputParameters.dat")
    def setInput(self, inputName):
        self.inputScript = inputName

    def setExecDir(self, execDir):
        self.exeDir = execDir

    # Organizes the command line to match the nucleation code's input style
    # ${execDir}/generalNucleation ${ScriptName}
    def organizeInput(self):
        self.clearProgArg()
        self.addProgramArg(self.exeDir+"generalNucleation")
        self.addProgramArg(self.inputScript)

#==================================================================================
class Towhee_Sim(Simulation):
    def __init__(self):
        Simulation.__init__(self, name="Towhee")
        self.setWorkDir(os.getcwd())
        self.addProgramArg("towhee")




#==================================================================================
if __name__ == "__main__":
    path = os.getcwd()
    print(path)

    nucTest = Nucleation_Sim()
    stderr = nucTest.runProgram_MPI(echoOutput=True, nProc=3)
#    stderr = nucTest.runProgram(echoOutput=True)

#    lammpsTest = LAMMPS_Sim()
#    lammpsTest.setWorkDir("/home/troy/SimulationCode/Tersoff_LAMMPS")
#    lammpsTest.setInput("/home/troy/SimulationCode/Tersoff_LAMMPS/in.freeze")
#    stderr = lammpsTest.runProgram_OMP(echoOutput=True, nProc=3)

#    towheeTest = Towhee_Sim()
#    stderr = towheeTest.runProgram(echoOutput=True)
#==================================================================================
