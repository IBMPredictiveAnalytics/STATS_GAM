* Encoding: UTF-8.
* basic tests.

get file="c:/spss24/samples/english/employee data.sav".
dataset name data.


stats gam dependent=salary variables=jobtime id=id errordist=gaussian
/nonlinear1 smoother=s variables =  educ prevexp parm1=6
/nonlinear3 smoother= bs variables = bdate parm1=6 parm2=3.

* error: missing parameter.
stats gam dependent=salary errordist=gamma
/nonlinear1 smoother=s variables=prevexp parm1=10
/nonlinear2 smoother=lo variables=jobtime parm1=.5.

* domain error.
stats gam dependent=salary errordist=gamma
/nonlinear1 smoother=s variables=prevexp parm1=10
/nonlinear2 smoother=lo variables=jobtime parm1=.5 parm2=8.

stats gam dependent=salary errordist=gamma
/nonlinear1 smoother=s variables=prevexp parm1=2.

