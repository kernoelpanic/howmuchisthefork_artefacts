--- mdp_double_spend.py	2022-03-09 12:54:55.008004162 +0100
+++ mdp_double_spend_modified.py	2022-03-09 12:55:32.759825543 +0100
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
