source_dir("datamakers")


#Now, for each scenario create an element of scenarios of the following form
add_scenario(dsc_znull,name="znull.1000",
            fn=null_z_datamaker,
            args=list(
              nsamp=1000
            ),
            seed=1:500)

add_scenario(dsc_znull,name="znull.100",
            fn=null_z_datamaker,
            args=list(
              nsamp=100
            ),
            seed=1:500)
