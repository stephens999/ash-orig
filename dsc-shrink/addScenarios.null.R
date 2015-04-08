sourceDir("datamakers")


#Now, for each scenario create an element of scenarios of the following form
addScenario(dsc_znull,name="znull.1000",
            fn=null_z_datamaker,
            args=list(
              nsamp=1000
            ),
            seed=1:500)

addScenario(dsc_znull,name="znull.100",
            fn=null_z_datamaker,
            args=list(
              nsamp=100
            ),
            seed=1:500)