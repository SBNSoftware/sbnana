cet_make_exec( testEarthModel
               SOURCE testEarthModel.cc
               LIBRARIES OscLibFunc
	     )
cet_make_exec( testOsc
               SOURCE testOsc.cc
               LIBRARIES OscLibFunc
	     )
cet_make_exec( testOscRndm
               SOURCE testOscRndm.cc
               LIBRARIES OscLibFunc
	     )

install_headers()
install_source()
