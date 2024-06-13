#!/usr/bin/env python
from numpy import *
data_1 = [
array([34.29885587793553]),
array([29.474738875120742, 30.396013428230525, 32.1179893134323, 31.520228661742763, 28.930037526336196, 30.30878721615335, 31.543776269350413, 31.030374352547994, 32.17919319385428, 29.370910119495324]),
array([32.25986660631723, 31.9475379476467]),
array([31.804102291312738, 31.963150278048097]),
array([28.42767451911602, 30.850265621479185, 25.99088103364648]),
array([30.115903997689003, 30.85747155424005, 30.923430815276404]),
array([28.250139843168107, 28.224021228661982, 27.744161678407174, 25.099086215354472, 24.42607977988648, 23.823829903778044, 24.732315275687725, 24.7087413810174, 25.188434118362753, 25.5212128808377, 23.67806880900789, 26.348309466556334, 23.533059521169655, 24.301510404023126, 27.092851947227423, 26.1915739452025, 27.930478328501646, 25.096182521142385, 23.419407578955457, 24.944361200411766, 24.821783589420775, 24.16517595510938, 27.26703950383811, 28.096629821872202, 25.541746194741084, 26.488940729597424, 27.03140824026795, 26.0321739549288, 24.414961170806386, 23.606873958838065, 25.196189308549357, 25.233345487616774, 23.951378752134236, 24.02791461865929, 23.88646229595341, 26.2430463799202, 25.733095936318474, 26.585256241128075, 23.146043500063506, 28.032495606936134, 25.72562074308543, 23.482630547870844, 23.695025075405994, 25.16317421708333, 27.741965176098834, 25.957372749148973, 25.398060548307317, 27.240474389082532, 25.96895914537459]),
array([27.825344650236755, 26.19792608761143, 25.966394693668654, 27.702129933965363]),
array([27.547374421532567]),
array([27.132406737221846, 26.842448336354874]),
array([23.79806876248672, 24.335226678886553, 23.377404154948927, 21.648064257196623, 20.323311013948143, 20.05422485804695, 20.423059708671165]),
array([26.024886737205353]),
array([24.725823272376484, 23.444963006096017, 23.472232767445753, 24.297746810153853, 25.2893093604401, 23.396745963165067, 22.145625538021527]),
array([24.197963010986708, 23.72693893360185, 23.186164112776407, 23.09559713275876, 22.553063956343713, 21.41366083944049, 22.013039334620743, 21.57395476012494, 19.553954913915035, 16.117135050633294, 19.560466078354363, 19.55709679731592, 18.343160311635092, 17.47027932964735, 16.669581026116127, 16.467983555849703, 13.91643123069132, 15.263374624641429, 15.589886506342703, 15.253933844251831, 14.717805126162075, 14.33192673444534, 15.34211242480878, 15.684256386671784, 14.725150384773979, 15.515635703934453, 11.456066562807996, 11.013808704701207, 7.4775651898127125, 10.32228542062856, 8.567384130108355, 8.527385434040406, 6.4660132171386495, 6.236083735883447, 5.380574307827175, 5.804173905509828, 4.047585777616153, 3.774744039806241, 0.0]),
array([24.312880689840064]),
array([23.59922919432982, 24.281988869158884, 24.40746402931533, 24.076750048375263, 23.128196413688737, 24.065981820536763, 21.91549880494994, 20.73562605917328, 22.79349466958386, 20.86082414647509, 21.876580620008788, 22.302650956034615, 21.484048358381862, 22.110148430382413, 21.477879228820235, 21.287894637537313, 22.78610866354261, 21.051716363264035, 22.831307429883086, 21.024370355068584, 22.091374869453958]),
array([21.978571047834524, 21.9647275616685, 20.30893622091472, 19.592365132789435, 19.466912190138864, 20.189281743626136, 17.55542437350874, 18.494818526917395, 16.063366957271555, 18.92484039877047, 20.34636659461247, 16.81817580942947, 14.998641499836074, 14.414922643222956, 15.659783458550011, 13.87217321123987, 14.081097389299515, 11.64450882409538, 8.032060253588936, 9.820208152320532, 8.191057561893071, 6.829912182049815, 7.228332376353367, 5.7679224312488575, 5.970276692001696, 6.128465023556654]),
array([22.326333829959534, 22.858704912135273, 21.960223867466254, 20.70986561065729, 16.959288454657624, 19.482502154112577, 18.428084283322278, 18.742891036871494]),
array([21.740663703119637, 21.92698697110477, 21.000818887845153, 21.3554219340395, 21.052283117717526, 21.851108186518136, 22.08173157299246]),
array([21.181771146989867, 20.713113076480628, 20.623827884822944, 21.60602831926176, 21.395984280436508, 21.33951712383922, 20.387816227798556, 18.706614387004574, 19.739581948814184, 19.547792371585174]),
array([21.55742493028014, 21.065536161286026, 20.645320195026777, 20.979713787516328, 20.978425739656178, 20.913172902995928, 21.55589974375965, 17.4814639569791, 19.603587238962465, 19.235939047724234, 18.32118336167247, 19.315246707863288, 17.444595074050618, 19.198412586044945, 19.901635342716276, 17.511873656521452, 18.42433195417445, 19.362338269412557, 17.349231437781125, 17.330387982400424, 18.438861288452692, 18.017311823856268, 18.73673759396913, 18.19238087194329, 17.481722001600943, 18.30840865807922, 19.716286132488655, 18.16613533785423, 19.21512735293675, 19.991450935115918, 20.235637237033064, 17.946230635195274, 17.28884941771726, 17.672207883799125]),
array([20.883358593585672, 21.011539687597672, 19.162348226405225, 19.240911304462013, 20.03215604393144, 19.26942885182715, 19.599843307465072]),
array([18.417430115512374, 19.121880094089846, 19.065566448148306, 18.054599630879856, 16.01783133272021, 19.943028562778558, 20.0754413616768, 16.438493663535702, 16.871051384410656, 17.485385364621933, 15.957121794077684, 14.000552905874873, 15.5334457234104, 15.83363047164935, 13.965424401819899]),
array([20.136142852514944, 19.358302522945703, 19.11667494375745]),
array([18.414086396929267, 16.61602388600689, 16.609642413365, 19.3577087946218, 17.693649631214804, 17.464180153873844, 18.40249860701952, 19.50584431742409, 19.153493703969723, 19.414683523578123, 15.997118403516309, 18.670469692388203, 15.274551648596981, 14.504583934742755, 14.598213146913952, 15.802076829919312]),
array([18.003661418050363, 19.61046995438232, 18.912110101553136, 19.177290251248525, 16.36004122714425, 18.34672550718987, 16.752095524924492, 19.804430818314852, 18.77838749866995, 17.51340774726694, 18.538799606785393, 14.324540205710084, 15.865399873327133, 14.663724223069778, 14.199501620618477, 15.688519377528293, 15.204765941535435, 15.099529975258998, 12.524937829797919, 13.431788255670307, 12.712340859851594, 12.674315959708093, 8.264274920971948, 10.578650870672455, 10.47399006213794, 8.197051085373856, 8.790317004737885, 9.362926512877062, 8.81508565876816, 11.419862061670585, 9.16840900242945]),
array([19.992770797486827, 19.69951513314455, 19.71355746640535]),
array([18.622579352156013, 18.826865797910237, 17.40616388637762, 17.732586279514933, 17.270858785095232, 19.313638921173386, 16.409723089756685, 17.654737031853625, 17.911913560427795, 18.379992350923587, 17.812087978780777, 16.402879603530376, 16.746240547102776, 16.91948578578802, 16.927191479593976, 17.418412112216053, 17.532073833292056, 18.010435380688527, 16.30987096056741, 16.089527191040105, 16.72167964438689, 17.203556284222373, 16.9631815237514, 16.344610208798894, 16.211687255165078, 17.863350182098856, 18.190870605651888, 19.472159924710695, 18.292575675566603, 17.110174344642033, 17.773777206206116, 18.772790288841062, 19.506310309222087, 15.7531079557533, 15.609232794425859, 15.706522627997623, 15.888213108373328]),
array([17.26645195760165, 18.11837799491786, 18.102273711284777, 17.62930840898292, 18.07701168253581, 18.941191865611053, 19.317132683339217, 16.42332249266788, 19.331568897117506]),
array([18.249891802520015, 18.857979230309958, 17.953244487086344, 17.993418692133638, 19.150993812323005, 18.196827692324934, 18.593726420095503, 18.567104215527166]),
array([17.803553408606437]),
array([16.88189649607708, 17.990338317873842, 17.13162315555766, 17.89235355525453]),
array([18.121912969737355, 17.482889965530592, 17.606429537201805, 16.686463525259434]),
array([17.770437498974122, 17.3365277814876, 17.546996455135503, 18.130080403041333, 17.74480232232571, 17.19476626235628, 18.02704604818925, 17.862572835126233, 17.30642809691682, 16.94617657575657, 17.84075940236409, 17.75168783467221, 17.870398245134965]),
array([18.16035888237727, 16.937226478973596, 17.22864816784405]),
array([16.979529231280225, 17.19923336541374, 16.198795100276545]),
array([17.32767938367074, 16.50189896709623, 16.060595451175892, 16.41459970180962, 17.08545257437763, 16.243280877939096, 18.11082619400218, 17.85477790493367, 18.05467946271799, 16.048099170304422, 17.228652848919527, 16.407970080895502, 16.285624259837782, 16.166212206984557, 17.70797130347943, 16.02396525817883, 16.78362500983978, 17.358720998833217, 14.55097584237124, 15.90094920026102, 15.380931060519984, 15.92995613049166, 14.187573609618807, 14.845893002503173, 14.440400423576065, 15.370907100295595, 14.669954296524226, 15.02009052964442, 14.642680319095643, 13.855361489955936, 15.883575288059458, 14.557430415314911, 15.013595962951916, 14.026529623663725, 15.142915143366492, 15.530478029281907, 15.785597378772732, 12.726494003236974, 13.740362888764128, 13.666685830626589, 12.887015685338705, 13.340391008703486, 10.866554989211464, 10.547866383523914, 10.76354560922457, 10.517038930144853]),
array([16.669518767973198]),
array([17.729085390681472, 16.801892337073685, 16.151507172370092, 16.390727866427813, 17.94794354363631, 17.101837080060275, 15.021744544560944, 15.494119808891334, 15.76749117902984, 15.077107681602115, 14.483836967171932, 14.979445251909334, 12.580163311147246, 13.40651645355157, 13.646413842516308]),
array([16.94329403470821, 16.197315684644263, 16.063870692971136, 16.335021905211708, 16.1000667935957, 17.474541968269502, 17.006865538340932, 17.132495972988433, 16.27361973713894, 15.617479906701996, 15.939231347007118, 14.547088530325555, 14.751141077487837, 14.426471613216744, 15.83002886170193]),
array([16.440544813693833, 16.25016490810797]),
array([16.02251238456398, 16.331932349831167, 16.72073921298254, 16.629821451791003, 14.642910333792015, 15.527338062355406, 13.89396486944823, 14.521387491733558, 15.742729434774823, 15.444196271467463, 14.739378484256639, 12.445142715272263, 13.673403423454243]),
array([16.23985551083852, 16.41251090945176, 16.133028330054408, 16.572950518406333, 15.977921285509685, 14.577499098438858, 14.991757314534006, 15.709467608001521, 13.957707729727387, 13.843969126835118, 14.466150259223005, 13.873682672796647, 15.736468389345692]),
array([14.504035887928607, 15.371602543535609, 14.29273587742988, 15.7330776638023, 15.319065146187679, 15.890073216750997, 13.525852929852215, 13.569934103943867, 9.521838977842956, 11.281408891608514]),
array([15.812056953697326]),
array([16.001335823426142, 15.354929347832552, 15.220471554134663, 14.13705634529109, 14.55255844784199, 15.779085059349006, 13.495343581551024, 13.81643734233313, 8.820598320031364, 9.583599314338185, 11.440545262922553, 7.91172404311866, 5.9047021557378905, 6.977879676861769, 5.542091221915147, 5.082605858395254, 3.9610802858616934, 4.872828858082019, 4.398955439554452, 4.472518428435532, 3.177847079530361, 2.250127534700075, 0.3588550086372445, 0.0]),
array([14.618429450178228, 15.575691191858144, 13.885881098023994, 15.564899069945765, 15.87884208938764, 14.177573335939755, 14.78339357043907, 13.853438450990708, 14.942071053057797, 15.45125815287335, 14.841466648376526, 14.541979993276035, 15.780104849851895, 13.16480359798539, 13.351104385913887, 13.58117505665967, 11.14207137032811, 10.277617530324534, 10.412783885564359]),
array([15.541246294127228, 13.97541788890615, 15.068475240163144]),
array([13.984679925033452, 14.690142321039554, 14.300231612543545, 15.378103807176045, 14.139482328613399, 14.05437707440479, 8.84442963832106, 10.328762865196062, 8.027673495088536, 11.314405302332496, 7.883601575292512, 9.933159163346158, 9.596375739011934, 6.179574175713707, 7.062067758965737, 6.876260079631148, 7.242756815286223, 5.625495034919148, 6.7156238914312585, 4.774255792637235, 4.176123112716024, 5.284682086953351, 4.254912069116137, 3.5629038840354723, 1.3652678706876453, 1.1693927720632102, 0.0]),
array([15.138689108623097, 15.11266732502282]),
array([14.702215641213561, 14.63061364130156, 14.579174354387538, 14.703763797365086, 14.790396374773863, 14.784359748880917]),
array([14.542290394438297, 14.311601738987525, 13.835179424121167, 14.099244365668628, 14.250304229777761, 14.348572386294299, 14.584246786919671, 13.498688811754597, 11.836989680306772, 13.432080892718336, 9.67574212796524, 7.601243295726484, 8.122794708398104, 11.293994425068139, 5.857860098154367, 6.824959265143439, 5.951863281102773, 6.34264665361648, 5.92098906708428, 6.576883618957122, 6.723994837036901, 5.4891280902456, 6.299769806273352, 3.7989558132351164, 3.2050840471242683, 2.476906996466916, 1.9260449024229331, 1.8681363317633015, 2.396049359551263, 0.6643285733054312, 0.4454889417685605, 0.023412768971273742, 0.0]),
array([13.932195630016569, 12.249473191693212, 13.204563132973362, 11.931603586155225, 8.925174288120013, 9.29544172209787, 10.786194333232931, 8.323879423293313, 8.668547366083727, 9.305027867562671, 11.16239911361649, 7.620715580820717, 6.747973096948012, 5.908997308910729, 5.418198277105572, 6.255470424233446, 6.716911376698812, 5.5893569992981895, 4.981835483344883, 5.143422974393576]),
array([14.131593800256942, 14.340837864026131, 14.106443888437862, 11.68562497364645, 10.387926560634993, 10.184472219068551, 8.071740787962588, 10.505024687530799]),
array([10.492107821851757, 10.385093078342567]),
array([12.042751975995824, 12.055940441594755, 7.389619158607217, 11.54687520446871, 7.280927918669421, 10.261070759664252, 11.080693804256542, 10.587965097030201, 8.041406644843594, 7.055703182696952, 6.241643520560081, 6.773726486126828, 5.54085902054533, 6.35973870142204, 6.841661239378516, 6.373018532015819, 6.623638278017375, 3.8104872102375564, 4.770671623121241, 4.72566319470567, 3.179618873455973, 2.216101598684055, 1.648660951565638, 1.2545456475110162, 0.24901084368346527, 0.07450846319985044, 0.0]),
array([10.431889315646131, 7.755541613607244, 5.397557390607259, 6.709507695681064, 5.573445940315757, 0.0]),
array([11.80631944979261, 7.970250444420207, 9.45788763806964, 11.220234923873214, 8.771789444754994, 10.852173481432493, 8.9982355542633, 9.453268618674999, 7.204101622031368, 6.154902109197004, 6.60445710029236, 6.576442975954642, 6.842905093140217, 5.672624511009798, 4.573278707085441]),
array([12.072045573622905, 11.793607220462057, 11.379050969609116, 11.468684605627832, 9.273396998733947, 9.198631755337825, 10.03340569773476, 9.398521955119113, 10.50303003250264, 11.439097494244626, 10.4500222544976, 11.047509186562694]),
array([12.12680379663202, 12.125413662493143]),
array([10.71250673091972, 2.481477940951338, 0.0]),
array([8.94945273156881, 7.267135299506309, 8.601943668443425, 9.573886808406348, 6.719381450399175]),
array([8.470702881708576, 8.74168551474584, 5.708625637477402, 6.740698245238363, 4.549175759193836, 2.7175346741133075, 2.996321422486408, 2.9895412629192837, 2.1299691705785153, 2.2266538584814466, 1.715189775271409, 1.1511884195940008, 0.1888576210757298, 0.3544766615162098, 0.09981345195472716, 0.0]),
array([9.003856948166167, 10.959475475279763, 9.74737381512405, 10.813736935943828, 9.446491750875868, 10.411585270666095, 10.546167387880185, 9.493454648914401]),
array([7.06856806868652, 6.6335665161257875, 6.32537973219202, 3.6405574663749096, 3.7108274268963637, 0.0]),
array([6.762297946556503, 0.34489897257291474, 0.0]),
array([8.055204439021425, 9.738243883147867, 7.697058823323699, 8.365632877250277, 10.035868555804525, 8.70383021365209, 7.428617367118385, 8.352195354596603, 9.886946600092518, 5.484961412101606, 6.320793868796434, 5.342859673839952, 5.476798832485839, 5.0179681713535205]),
array([10.279067989505286, 9.763723869568775, 9.415737140700024]),
array([8.652822151934702]),
array([8.037136403524876, 9.965079340497393, 7.065514428075355, 5.84082051617961, 6.554812193551172, 7.192490206146916, 6.1313928041004, 5.755273593590769]),
array([10.210700768431822, 7.87657771495717, 8.783837965914797, 9.593546932557258]),
array([9.578301245179997, 9.547927549000041, 8.59783244385643, 9.939373443025662, 9.063484880657263, 9.329190196703804, 6.835500297180741, 6.40286281818886, 6.320321515440035, 4.437345860482914, 2.653496218631243, 3.407297459925869, 3.133682201520028, 2.090931010463211, 2.347883852040348, 0.6868682146120682, 0.0759914512637956, 0.0]),
array([9.275320429123994, 7.681524137842917, 9.084110614229894, 9.634280901481144, 8.59233708113095, 7.411190975477367, 7.9671976892569685, 9.664306179723802, 7.926011233709751, 7.811428707995663, 9.648777105029401, 8.344958160682104]),
array([9.185165392989157, 8.468552168248587, 9.320288618650858, 6.943844954506852, 6.696032925103275, 4.788840652747137, 3.2888683124746003, 3.111173090378866, 1.8921640436775655, 0.05978685449641363, 0.0]),
array([8.981228508448353, 8.49649861908309, 8.345776947748435, 7.8981708983983, 7.742490762703594, 7.014247031082978, 6.328133960871704, 6.815213865555046, 5.932879722147995, 5.854581112163462, 6.772906900383874, 5.131743975717048]),
array([8.120976883722632, 9.257289056298777, 7.969315673149336, 8.500307402770712]),
array([8.644609931106663, 9.159532286522218, 9.240940476510007]),
array([7.902113748022141, 7.543240390463127, 6.348170587638556, 5.982475328068709, 7.2338217752953105, 5.879539791709414, 4.397645871347702, 4.250473605126271, 3.210097500060761, 3.507049411327931, 2.3253537402141555, 2.4528482466166475, 1.7649007671758097]),
array([7.5064496255908395, 8.570165782296886, 7.679858253929854, 7.8268679894618085, 8.236734801815702, 7.458370588807485, 8.508614367844483, 6.992170918015664, 7.212341801765497, 6.843274989863777, 5.653342095391254, 7.0515495508854915, 5.766890445983998, 7.112563285879992, 5.7807364119306275, 6.230493684649752, 6.889170340476667, 6.700630805223589, 6.76093954138788, 5.901487212131795, 7.215362838803048]),
array([8.331656380993758, 7.474716878391026, 8.358962580590415, 7.128029606521971, 6.505731833418018, 5.869692898061205, 6.2152194273113075, 5.923640050590768, 6.4101005019510175, 6.066605591624195, 6.61795714850733, 6.242840938115374, 5.3579747506115565, 7.190966357680115, 6.873727605469529, 6.006655946423759, 6.215574637635451, 6.183665003439593, 6.212844080460206, 6.28810297349927, 4.150251212188392, 3.822467788842178, 5.17424627057334, 4.367397701644622, 5.078841768487468, 4.717099621470591, 3.987329689015176, 4.892998878008617, 2.680749249650789, 3.137098966092929, 2.08523536537837, 2.2111886537587795, 2.3053307302911175, 1.7852030120546118, 1.020038711241503, 1.3461478535823659, 1.0544116721183152, 0.4584088497017979, 0.06984891029919525, 0.03621922519944605, 0.0]),
array([8.404944772488811, 7.778735605656266, 7.659468783745865]),
array([6.778770922695362]),
array([8.05783965700743, 8.368875414702654, 7.990952344428564, 7.819735918021495, 8.286409382524955, 6.5698115143698095, 6.988117127171251, 7.128053883564905, 6.948143312980782, 7.015536612187841, 6.72869920242239, 6.391580770602184, 6.959951060374037, 6.680890273544687, 6.67215774607189, 7.076968130906426, 6.4543144817031335, 6.897428472414171, 6.910559993700766, 6.545954553914771, 6.5671848164870985, 6.923960527347008]),
array([5.589409560639808, 6.070592596543396, 6.29863293905065, 5.427792686805367]),
array([6.981830967383247, 6.957202101155935, 7.0548799800218465, 6.121773372738458, 6.230742965042257]),
array([8.14765277769244, 6.348358615352784, 5.872759327804521, 6.351343011290203, 6.061855072273212, 6.947612030707728]),
array([7.665140653566586, 7.5203235535571755]),
array([7.502078099552986]),
array([8.051590195327448, 6.934364776336594]),
array([7.734145250088609]),
array([7.877767588997076, 7.833921702318247, 7.438030410645952, 7.7263032731031185, 5.578633364122469, 6.817260977046301, 6.236000145419362, 4.893964058014867, 5.0122470382389945, 2.6578420699642162, 2.302038268115742, 2.465735974341591, 2.415932505404135]),
array([7.017279960849455]),
array([7.480372637988937, 7.504819907370709, 7.8572196980899385, 6.50070976857286, 6.961509673897431, 5.382547833630641, 7.105205138274203, 6.675600488638454, 7.173719668071294, 5.5403805855000785, 3.9413395037091936, 4.76953415015193, 3.6410513616465674, 3.451008494516455, 2.1173341940248536, 1.9977298316337415, 2.434630637064574, 2.1353108942946624]),
array([7.418467121085168, 7.868410192851755, 6.628676991497324, 5.4451108316969234, 5.568400235822478, 6.66556280272874, 6.5320074309030085, 5.957508711523804, 6.509925689023334, 5.977078313011514, 7.143907732061533, 4.571743177447337, 3.5707627168797837, 3.123438853991869, 3.3316968953833648, 2.157095440498878, 2.5532578466744957, 2.3148078397432417, 2.3018960689046435, 1.1264072564986352, 0.32617130120220467, 0.08055264932432704, 0.06583217564419445, 0.0862038492766686, 0.0]),
array([5.5587883077642095, 5.107971059055783, 0.9803804141223459, 0.0]),
array([5.824270570639697, 5.577599608958032]),
array([6.742543584828952]),
array([7.667567737131692]),
array([6.106892590436804]),
array([6.950394479265608, 7.089176615092588]),
array([5.431170526388552, 6.295603645198284, 6.916022633182216, 6.572409716798435, 4.038771338173279]),
array([6.228795749908385, 5.866782161774029, 5.550526087433532, 7.063055279983165, 7.081905390914548, 7.075712860811399, 4.65290066016376]),
array([6.2714932193344755, 5.638953885613832, 6.430642410298511, 6.297867812558908, 5.4426460573203475, 6.369755828770337, 6.704551784145388, 6.3572045306393274, 5.8763093560291075, 4.900255494768893, 2.02458517721617, 0.4249418137003572, 0.1938359871153036, 0.0]),
array([6.6687520638531534, 6.725126029673755, 1.9761191376183023, 2.185163792572482]),
array([6.620651005386953, 6.061806577335201, 2.7827846564894188, 3.1674628844357415, 2.0477694704431038, 0.6203916407842656, 0.36657134225882904, 0.03182861737413162, 0.0]),
array([6.1996997674004986, 5.667958799522718, 5.623537908036709, 5.46159194362313, 5.694345952569947, 6.221659230483682, 5.7063881595326915, 6.430187319270738, 5.90836401835549, 3.8771081378216365, 3.9826243932361667, 3.981030367361315, 3.4859424910291823]),
array([6.096507953659992]),
array([5.98144048851774, 5.498973543369308, 3.217713211133451, 2.575491641471194]),
array([5.387705019175346, 5.4083534759329135, 5.7575066133436446, 5.418182825357526, 4.73677658051159]),
array([5.818248295486424]),
array([5.339065245356011, 5.721563954191206, 5.434747246097213, 5.095847563433568, 4.911342888322121, 3.7178161207445486, 4.630275211831698, 2.945325192926926, 3.4559309015429776]),
array([5.763480830813274, 4.806446613979477, 4.0650647671249]),
array([5.34074919357483, 4.9474490669347375]),
array([5.483869743788569]),
array([4.140105888534391]),
array([2.458691615609759, 0.6095188262635005, 0.0]),
array([5.373220586407456, 5.0411917858366895, 2.1441266674156334]),
array([4.727887589602158, 4.223229075584284, 3.4384718197740267, 2.5028142257050368, 2.429891534172213, 2.548679667887113, 0.9900951301534728, 0.33688240386278645, 0.5775808015610328, 0.5694267613409842, 0.5696568651718641, 0.42121406496633546, 0.05328729894181264, 0.0]),
array([4.281133136672699, 5.165095307520767, 3.9250275753027477, 5.118535282280202, 3.5831731796919373, 2.1271677499971373, 1.269504242159564, 0.5828851824957946, 0.3915540110358006, 0.16935903716034073, 0.0]),
array([5.176847077734104, 3.8308738600915193]),
array([4.176302489333325]),
array([3.97618393699715, 4.168564125758554, 5.001413998001661, 4.699010416834017, 2.724874494099449, 2.260140903022142]),
array([2.6926591315201702, 2.575561458224556, 1.9873103988891263, 2.355376046892092, 1.1961225614912485, 0.6569188467420295, 0.0]),
array([4.273194352654319]),
array([2.396234675070649]),
array([4.71714144064965]),
array([4.30128484643327]),
array([4.1151275187724385, 4.180874286351596, 2.750062164769387, 2.4956561498866177, 2.5228445208440857, 2.1262782210443816]),
array([4.262565157497864, 3.8182815274274167, 1.8487404656358524, 2.538063798782524, 0.0]),
array([4.115491710860112, 4.027003902894941, 1.8276926842584413, 0.3418281656316034, 0.0763022194377513, 0.0]),
array([3.9284708379179336, 3.9068607005191085, 3.7625232438146385, 2.777406371367773, 2.663674093120873, 3.311443479604453, 3.4932134783021094, 3.5776100103524815, 2.4405028183564905]),
array([4.053108785112066, 4.021120882516029, 2.767911566307873, 2.3188963647880976, 2.030862736171916, 1.1734694084422372, 0.7488227509051346, 0.0]),
array([2.0914372111288366, 0.0]),
array([3.559093612553389, 2.849287372482941, 2.8638181353200216, 2.2127831962862765]),
array([3.763757785687339, 0.6894914297672687, 0.0]),
array([3.604773342734179, 2.897774158449994, 2.9003625001991167, 2.7083890092385348, 3.093930432487491, 3.355535947729047, 2.6077472180417507, 2.758264490854305, 2.456580582530622, 1.9617963828050848, 2.3702552169458775, 2.0579136207864868, 2.0887377625034786, 1.1545780599334616, 0.9031100574330428]),
array([2.142658113534875, 2.5624902195233705, 0.5525290936807222, 0.0]),
array([3.0617781744783157, 3.3630070661082354, 2.1310938310838305, 2.53816137117614, 0.1110427164292567, 0.0]),
array([2.7222940638736977, 2.464806045892009, 2.217775332477935, 0.7094756850647483]),
array([3.7972117999345176, 3.4861234089296276, 1.546310801083394, 0.0]),
array([1.9106700341904537, 2.186845582796949, 2.506522429774029, 1.5049228687589409, 0.4223129133400037, 0.7799476690672542, 0.0]),
array([2.826742073831367, 1.023223982336884, 1.394421189976859, 0.5907271550337434, 0.16202999193288747, 0.0]),
array([2.7440170754578252, 2.3901967347940416, 2.2905110099420622, 1.1381854295671001, 0.0]),
array([3.2178179814357115, 2.123100796019372, 2.379672588372961, 0.5270954709120443, 0.3312928449689715, 0.6254229980252523, 0.016148641341771977, 0.0]),
array([2.896247729235229, 2.001913392071579, 2.4999602501425646, 0.9428293942866935]),
array([2.228469978919346, 2.120984245671919, 1.1115484497576986, 0.0]),
array([2.3971041960444546, 2.343623818417468]),
array([1.9819970352044987, 2.2108127504678308, 0.0]),
array([2.444883370822524]),
array([2.194571552659512, 2.489787561000832, 1.845987086872845, 2.0052206461932407, 0.9876549269752944, 0.8996426881278695, 0.8905505195730525, 1.3224754181528433, 0.5996943529015535, 0.3418656311540575]),
array([1.816702384726053, 1.9258656684201314, 0.7655588584818338]),
array([2.266619194890172, 0.06220735292806849, 0.0]),
array([2.265507676677819, 2.055561384499536, 1.2751391365822995, 0.22407133710756988, 0.07969970956189446, 0.0]),
array([2.7938935431826866, 2.806903521369009]),
array([2.244352300988245, 2.5015037079615126]),
array([0.3461658105334825, 0.47349213352781394, 0.0]),
array([1.9495926616326886, 0.0]),
array([2.0271145722619286, 2.3486523165033035, 0.0]),
array([1.974564306990337, 2.0265227827763956, 0.0]),
array([2.4614793342521835]),
array([2.26302749302526, 0.0]),
array([0.591091112902667]),
array([2.053143637832581, 0.5483794947359248, 0.0389756174661606, 0.0]),
array([2.12973477796127]),
array([0.523603200528458, 0.5042470332969238, 0.005217480130141741, 0.04767339227657623, 0.0]),
array([1.4524924216858741, 1.6859094165372646, 1.6034032249388979, 0.86614727985268, 0.47932114813011106, 0.10494662859646176, 0.0]),
array([1.5103970384745176, 0.06695148262397793, 0.035507816573675896, 0.0]),
array([0.9421900459387932, 0.09991638334714517, 0.0]),
array([1.3091812085635803, 1.5786106203471562, 0.246248883422998, 0.0]),
array([0.8443283530756417]),
array([1.6744678088036242, 0.10297956936079479, 0.0]),
array([0.5098177820679921, 0.0]),
array([1.3413043101417452, 0.5514101580166031, 0.0]),
array([0.12178931008118646, 0.03411862380642923, 0.0]),
array([0.2861267767391585, 0.0]),
array([0.24710865660866832, 0.0]),
array([1.0551071405484622, 0.44201875607405144]),
array([0.2104526041657906, 0.31509901582224104, 0.6153362071337062, 0.0]),
array([1.2043919925977793, 0.0]),
array([0.05617337685431316, 0.0]),
array([1.073052477619898]),
array([0.8426229463142165, 0.0]),
array([0.338411465077996, 0.0]),
array([0.9220732716599583, 0.0]),
array([1.0061519698950785, 0.1711879258485821, 0.08302208107673159, 0.031919245450349826, 0.0]),
array([0.1607205843805215, 0.05212336301567955, 0.07391998942349061, 0.0]),
array([0.9558976751153295]),
array([0.5732790029073983]),
array([0.31297623940654357, 0.0]),
array([0.10696008945647058, 0.0]),
array([0.1991333666643883, 0.0]),
array([0.6173207088395134, 0.0]),
array([0.10919750766199918, 0.0]),
array([0.1114029971040808, 0.0]),
array([0.026761573236734643, 0.0])
]
d = [data_1]
names = ["83"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T1', 'T3', 'T4', 'T6', 'T7', 'T8', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'T36', 'T37', 'T39', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T46', 'T47', 'T48', 'T49', 'T51', 'T52', 'T53', 'T54', 'T55', 'T56', 'T57', 'T58', 'T59', 'T61', 'T62', 'T63', 'T65', 'T66', 'T67', 'T69', 'T70', 'T71', 'T72', 'T74', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T85', 'T86', 'T87', 'T89', 'T90', 'T92', 'T93', 'T94', 'T95', 'T97', 'T98', 'T99', 'T100', 'T101', 'T102', 'T103', 'T104', 'T105', 'T107', 'T109', 'T110', 'T112', 'T113', 'T114', 'T115', 'T117', 'T118', 'T119', 'T120', 'T121', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T129', 'T131', 'T132', 'T133', 'T134', 'T135', 'T136', 'T137', 'T138', 'T141', 'T142', 'T144', 'T145', 'T147', 'T148', 'T149', 'T150', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'T158', 'T159', 'T162', 'T163', 'T164', 'T168', 'T169', 'T170', 'T171', 'T172', 'T173', 'T175', 'T178', 'T179', 'T180', 'T181', 'T182', 'T183', 'T185', 'T186', 'T187', 'T188', 'T190', 'T191', 'T192', 'T193', 'T195', 'T197', 'T199', 'T200', 'T202', 'T204', 'T205', 'T206', 'T207', 'T209', 'T210', 'T215', 'T216', 'T217', 'T218', 'T219', 'T220', 'T221', 'T222', 'T223', 'T225', 'T226', 'T227', 'T228', 'T230', 'T231', 'T232', 'T233', 'T234', 'T236', 'T239', 'T241', 'T242', 'T243', 'T244', 'T245', 'T250', 'T255', 'T257', 'T259']
def get_taxa_names(): return taxa_names