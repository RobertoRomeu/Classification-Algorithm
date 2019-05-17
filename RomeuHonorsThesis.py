import time
import os
import math
#import parselmouth
"""
tuginoaiha sigatu 8tu nijihanni aou

Draga-kara
    (1) Figure out how to port to Cornell box
    (4) Background-wo kaki-hajimeru
        What are similar types of things j fill
        that in and distinguish
    References earlier
    Email Draga draft for structure/explanation

    email mats rooth about being second reader
    soc net corp
    draft
"""


"""Static variables"""
mes =5      #Number of measures
pp = 10     #Number of vertical partitions on total data
hpp = 5     #Number of horizontal partitions

class spkr:
    def __init__(self, n,raw=True):
        if raw:
            self.id = n
            s="s"+str(n)+"form.txt"
            d = fyomu(numbPath(n),s)
            self.main = mOut(d)
            self.spl = split(d)
            self.parts =  giveParts(d)
        else:
            self.main = mOut(d)
            self.spl = split(d)

    def compTo(self,spkcomp):
        total = []
        for i in range(0,mes):
            total.append([])
            total[i].append(tval(self.main[i],spkcomp.main[i]))
            for j in range(0,pp):
                total[i].append(tval(self.spl[j][i],spkcomp.spl[j][i]))
        return total

class weigh:
    def __init__(self, pool):
        self.pool = pool
        self.ws = getMWs(pool)
        self.reject = 1000000.0      #Rejection mechanism to implement

    def wval(self,x,y):
        s=x.compTo(y)
        sum=0.0
        for k in range(0,mes):
            for l in range(0,pp+1):
                sum+=self.ws[k*(pp+1)+l]*abs(s[k][l])
        return sum

    def guess(self,sp):
        min=self.reject
        g=None
        for x in self.pool:
            y=self.wval(x,sp)
            if y<min:
                g=x
                min=y
        if min==self.reject:
            return "Not in Pool"
        return g.id

    def cmat(self,pool):
        cm = []
        ll=len(pool)
        for i in range(0,ll):
            cm.append([])
            for j in range(0,ll):
                cm[i].append(0.0)
        for i in range(0,ll):
            for j in pool[i].parts:
                res = self.guess(j)
                if res != "Not in Pool":
                    cm[i][int(res)-1]+=1.0
        return cm

class cd:
    #idiom to use as context
    """Context manager for changing the current working directory
    Use: "with cd(path):"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def giveParts(d):
    """Returns horizontal partitions of raw value /d/, manifest as spkr list"""
    PP = hpp
    parts = []
    for i in range(0,PP):
        parts.append(spkr(part(d,PP,i),False))
    return parts

def getMWs(sp):
    """Returns a (//pp//+1)*//mes// list of weights (sum to 1)
    in the order main>partitions>next main
    Maximize those most inter-speaker different
    """
    ll = len(sp)
    ws = []
    for i in range(0,(pp+1)*mes):
        ws.append(0.0)
    for i in range(0,ll):
        for j in range(i,ll):
            s = sp[i].compTo(sp[j])
            for k in range(0,mes):
                for l in range(0,pp+1):
                    ws[k*(pp+1)+l]+=abs(s[k][l])
    return norm(ws)

def norm(lx):
    """Norms a list /lx/ by summing then dividing. Assumes /lx/ is only
    summables
    Returns the normed list"""
    sum =0.0
    for x in lx:
        sum+=x
    ws=[]
    for i in range(0,len(lx)):
        ws.append(lx[i]/sum)
    return ws

def split(s):
    """Constructs the pp vertical partition of a data stream s,
    returns a list of each partition with the measures'
    list of formant [mu,stdev,n], levels to a 1-norm of the step"""
    max = []
    min = []
    for i in range(0,mes):
        max.append(-60000.0)
        min.append(60000.0)
    for j in range(0,mes):
        for x in s[1][j]:
            if(x != None):
                if x > max[j]:
                    max[j]=x
                elif x < min[j]:
                    min[j]=x
    #defines the max-min
    step = []
    for i in range(0,mes):
        step.append((max[i]-min[i])/float(pp))
    #defines the step for each mes
    end = []
    for i in range(0,pp):
        end.append([])
    for j in range(0,mes):
        parts = []
        weys = []
        for i in range(0,pp):
            parts.append([])
            weys.append([])
        for i in range(0,len(s[0])):
            if s[1][j][i] != None:
                k = chkAdd(min[j],step[j],s[1][j][i])
                parts[k].append((s[1][j][i]-min[j]-k*step[j])/step[j])
                weys[k].append(s[0][i])
        for i in range(0,pp):
            if len(weys[i])>1:
                end[i].append(wsamp(weys[i],parts[i]))
            else:
                end[i].append([None,None,0])
    return end

def chkAdd(min,step,x):
    """Helper for vertical step calc, takes the /min/ of the range,
    has the /step/ size and an observed /x/
    Returns index on //ws//"""
    for i in range(0,pp-1):
        if (x >= min + step*i) and (x < min + step*(i+1)):
            return i
    return pp-1

def fyomu(pathIn,n):
    """Reads the /n/-th speaker data of the processed WAV file
    from the Buckeye corpus. Specifically, reads for "To Burg" formant
    analysis
    Outputs a list with i-0 being the time-step intenstity, and i-1 is a list of
    the time-step for formants. If a formant is not present, it's value is None"""
    with cd(pathIn):
        txt = open(n,"r")
    txt = txt.read()
    si = 0
    for i in range(1,6):
        si = txt.index("\n",si+1)
    ni = int(txt[si+1:txt.index("\n",si+1)])
        #Count n of windows
    for i in range(5,9):
        si = txt.index("\n",si+1)
    strm = []
    form = [[], [], [], [], []]
    for i in range(0,ni):
        strm.append(float(txt[si+1:txt.index("\n",si+1)]))
        si = txt.index("\n",si+1)
        ct = int(txt[si+1:txt.index("\n",si+1)])
        si = txt.index("\n",si+1)
        for j in range (0,5):
            if j < ct:
                form[j].append(float(txt[si+1:txt.index("\n",si+1)]))
                si = txt.index("\n",si+1)
                si = txt.index("\n",si+1)
            else:
                form[j].append(None)
    return [strm,form]

def wsamp(w,sitta):
    """Takes a list /w/ of weights, index-linked to
    observed measures /sitta/ and returns the sample mean,
    st. dev, & size (n) of the observed"""
    n = len(w)
    emp = set()
    mu = 0.0
    wei = 0.0
    for i in range(0,n):
        if sitta[i]==None:
            emp.add(i)
        else:
            wei+=w[i]
            mu+=w[i]*sitta[i]
    mu = mu/wei
    m=n-len(emp)
    sd = 0.0
    for i in range(0,n):
        if i not in emp:
            sd+=w[i]*((sitta[i]-mu)*(sitta[i]-mu))
    sd = sd / ( ((m-1)/m)*wei )
    sd = math.sqrt(sd)
    return [mu,sd,m]

def tval(x1,x2):
    """Returns Welsh's t-test distance of two samples"""
    if x1[2]==0 or x2[2]==0:
        return 00.0      #Need to resolve later what to do when tiny
    t = x1[0]-x2[0]
    und = ((x1[1]*x1[1])/(x1[2])) + ((x2[1]*x2[1])/(x2[2]))
    und = math.sqrt(und)
    return t/und

def part(x,par,ind):
    """Partitions a raw /x/ into /par/ partitions, returns the
    /ind/-th partition"""
    n = len(x[0])
    ns = int(n/par)
    strm = x[0][ind*ns:(ind+1)*ns]
    sub =[]
    for i in range(0,mes):  #mes here
        sub.append(x[1][i][ind*ns:(ind+1)*ns])
    return [strm,sub]

def mOut(fyomi):
    """Returns the synopsis of wsamp (mu,stdev,nct) for a
    data set /fyomi/"""
    whole = []
    for i in range(0,mes):  #mes here
        whole.append(wsamp(fyomi[0],fyomi[1][i]))
    return whole

#Possibly not needed  or needs modified
def kakuComp(ref,compset,name):
    """Writes the t-distance to /ref/ from a sequence of
    /compset/ and writes it to the desktop as .txt file with
    name /name/

    """
    N = len(compset)
    s=""
    for i in range(0,mes):  #mes also here
        for j in range(0,N):
            s+=str(tval(ref[i],compset[j][i]))+"\t"
        s+="\n"
    file = open(name,"w")
    file.write(s)

def numbPath(n):
    """Returns the directory path to find the file
    (For Buckeye)"""
    s="C:\\Users\\rober\\Desktop\\CRUDO\\s" #Update path here
    if n<10:
        s+="0"+str(n)+"\\"
    else:
        s+=str(n)+"\\"
    return s

def odds(cm):
    """Takes a confusion matrix /cm/ and calculates
    the accuracy as proportion of guesses on the diagonal
    """
    ll=len(cm)
    sum=0.0
    for i in range(0,ll):
        sum+=cm[i][i]
    t=float(hpp*ll)
    return sum/t

def makeMeHappy():
    s="s"+str(3)+"form.txt"
    d = fyomu(numbPath(3),s)
    s= ["","","","","",""]
    for i in range(0,len(d[0])):
        if i%250==0:
            s[0]+=str(d[0][i])+"\t"
            s[1]+=str(d[1][0][i])+"\t"
            s[2]+=str(d[1][1][i])+"\t"
            s[3]+=str(d[1][2][i])+"\t"
            s[4]+=str(d[1][3][i])+"\t"
            s[5]+=str(d[1][4][i])+"\t"
    s=s[0]+"\n"+s[1]+"\n"+s[2]+"\n"+s[3]+"\n"+s[4]+"\n"+s[5]
    filef = open("250pt sample data.txt","w")
    filef.write(s)

def makeMeHappier():
    s="s"+str(3)+"form.txt"
    d = fyomu(numbPath(3),s)
    s = [ [], [] , [] , [] , [] , [] ]
    for t in s:
        for i in range(0,600):
            t.append(0)
    for x in d[0]:
        s[0][indexRange(x)]+=1
    for i in range(0,5):
        for x in d[1][i]:
            if x!=None:
                s[i+1][indexRange(x)]+=1
    st=""
    for i in range(0,600):
        for j in range(0,6):
            st+=str(s[j][i])+"\t"
        st+="\n"
    filef = open("histogram.txt","w")
    filef.write(st)

def hist2():
    s="s"+str(3)+"form.txt"
    d = fyomu(numbPath(3),s)
    s = []
    for i in range(0,50):
        s.append(0)
    for x in d[0]:
        s[iR(x)]+=1
    st=""
    for i in range(0,50):
        st+=str(s[i])+"\t"
    filef = open("histogram2.txt","w")
    filef.write(st)

def indexRange(x):
    return int(x/10.0)

def iR(x):
    return int(math.log(x)/math.log(.65))

"""Here we begin writing"""
print("\n\t\tWe start here")
#cmrate = []
#for i in [2,3,20]:
#    if False:
#        pp = i
#        cmrate.append(main(20))
#        print(time.process_time())
#makeMeHappier()
hist2()

#print (cmrate)
#pp=2
#main(20)
#pp=5
#main(20)
#pp=10
#main(20)
#Timer for process measurement
print ("\nTotal time: "+ str(time.process_time()))
