Notes on setting up a virtualenv
======================
Virtual environments for python are great on codon, because you don't have to bug IT every time you want to check out a cool new python module. But I always forget how I went about setting these up, even though they are very easy. 

I like to put my virtualenvs into a separate directory - This would make it easy if you ever decide you need a new one, or a specific set of modules you dont want to trample over. I keep my main one named 'base'

0. mkdir $HOME/virtualenvs/
1. module load python 
2. virtualenv --system-site-packages $HOME/virtualenvs/base

	the --system-site-packages lets you use the system versions of things if they are already installed. I believe you can use your own versions by explicitely installing them to your virtual env. 
3. source $HOME/virtualenvs/base/bin/activate # need to do this whenever you start up your term, can add it to your .bashrc if you like. 
	Don't forget to do this in your qsub scripts. I go module load python, source $HOME/virtualenvs/base/bin/activate

Now pip install -U whateveryouwant should work, and should install to your virtualenv. 

