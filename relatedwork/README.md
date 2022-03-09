# Double Spending MDP

This folder includes some related work on MDPs from Gervais et al. presented in https://eprint.iacr.org/2016/555.pdf
and forked from https://github.com/arthurgervais/pow_mdp 

To run the double spending MDP the steps to set up a python2 environment with the necessary packages are described for Ubuntu 20.04:

```
apt install python2 
apt install python-tk
pip uninstall virtualenv # If you get some error from virtualenv use the version in the apt repository
apt install virtualenv
virtualenv -p python2 venv2
source venv2/bin/activate
pip2 install -r requirements.txt
```

A small change has been made to describe and reflect the correct ordering of the command line parameters: 
```
$ diff -u mdp_double_spend.py mdp_double_spend_modified.py > mdp_double_spend_patch.py
...
@@ -443,12 +443,16 @@
         cost = float(sys.argv[1])
         gamma = float(sys.argv[2])
     else:
-        print "Not enough arguments"
-        print "Usage: %s <gamma> <cost>" %sys.argv[0]
+        print "Not enough arguments:"
+        print "https://eprint.iacr.org/2016/555.pdf"
+        print "gamma  = propagation ability, fraction of the network that receives advasary block first"
+        print "cost   = mining cost in terms of block reward, [0,\\alpha]"
+        print "Usage: %s <cost> <gamma>" %sys.argv[0]
         return
     stale = 0.0041
-    k = 6
-    hashrate_k_plot(stale, gamma, cost, cutoff=20):
+    #k = 6 # actually this variable has not effect ... 
+    #hashrate_k_plot(stale, gamma, cost, cutoff=20)
+    hashrate_k_plot(stale, gamma, cost, cutoff=5)
 
 if __name__=="__main__":
     main()
```
The modifed cutoff should yield faster results. 


Run for example with: 
```
$ time python2 mdp_double_spend_modified.py 1 0.15
```
