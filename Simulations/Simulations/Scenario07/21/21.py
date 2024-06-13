#!/usr/bin/env python
from numpy import *
data_1 = [
array([34.1057049977952, 34.45074308194396, 34.056587051316676, 34.464453129398336]),
array([32.329890564708634]),
array([31.30111197718249, 28.76920178015316, 31.397163430083705, 32.8329181689785, 28.474866506394925, 30.74437556011926]),
array([25.720938289402554, 25.555045646869278, 20.77759863042919]),
array([28.797780521458368, 28.03652933732185]),
array([24.14008328294433, 25.92002731853602, 23.314438879016944, 26.418511568948993, 25.311845910012973, 26.30948414858738, 22.122141611894087, 22.20714663488441, 22.59423528830473, 21.810294344279345, 21.9093647618454, 22.34906090732751, 21.849920842925926, 20.347921622021307]),
array([24.345805929735317, 24.891474550070996, 24.75340861104645, 24.978527542425518, 24.094950834882333, 24.57321068685287, 25.321316611320235, 22.337096860594535, 21.20752176578847, 22.73035319036825, 20.54325585231287, 22.385027949219744, 22.11993414495302, 21.561085244203184, 20.94300165292878, 22.50535229757866, 21.75154101856255, 21.44760182602574, 22.839520512054513, 21.709513308725413, 21.82156140512607, 18.00416188723382]),
array([23.031990457663593, 21.39897160418532, 20.734718666010785, 21.23974533046647, 20.502681189794863, 21.029286329170283, 22.01429942852375, 21.998279742686982, 21.839050072584257, 22.451852649765975, 22.556252533555906, 21.291824698140086, 21.995309239632007, 21.27927912277237, 22.94993078648494, 21.527672322578965, 22.06248177801369, 22.410155593396418, 22.34329536510534, 21.871746836797513, 22.534507080393624, 21.96278024909635, 18.376136999084736, 17.177895813084625, 17.499338904450045, 18.1493609769875, 15.413603542626795, 14.985089620312106, 15.155652557899513]),
array([20.664008266791548, 22.99165776053696, 22.744286257922845, 20.757360165636285, 20.122101618527253, 16.402253292468778, 17.303853092795915, 14.167876382736521, 11.97075161947867]),
array([21.545863724427456, 21.00474622096437, 20.458216921528493, 20.58289897619546, 21.386170689765144, 16.99425233351206, 19.714196432802996]),
array([21.195642447339797, 20.583699157134397, 20.9721757239273, 20.935857242211597, 21.3236214294062, 21.39975947715363, 19.31753928838281, 18.896962656451045, 19.050128818816262, 20.419007008464003, 15.907395575250183, 15.106357356995549, 14.224450489009154, 15.500692246676529]),
array([20.606762369210617]),
array([20.5363505024308, 17.753687569687713, 15.532575233213494, 14.440696479146618, 8.225340409283099, 8.133788016856899, 9.583248271069184, 8.724538678847294, 9.043670800642499, 6.212588769465993, 6.143517353530278, 5.811115563573626, 4.4588270892164426, 5.183916072425548, 3.2095607327221036, 3.163204716174358, 2.7538274008667805, 3.100098663635315, 1.5286852434357532, 1.2955958874926656, 1.410130567096432, 1.207754505592879, 0.0]),
array([20.48616490687536, 20.653657233969838]),
array([20.56336619105968, 16.7774176359112, 20.025427510287138, 14.192059783898308, 15.640786302264281, 14.04922319324849, 12.70906741853843, 8.911040205043072, 10.74253881498923, 11.376414867179102, 7.431115991729371, 9.813721624837145, 5.822971280038089, 5.820663460936903, 6.107192913727991, 5.764643823987191, 4.717756962946346, 3.3177071845361716, 3.067902468979619, 2.948931726998233, 2.5866322434464175, 3.595743013419085, 0.8130325226323344, 1.6440803536954025, 1.0552233237670787, 1.07777504648096, 1.7479860670020786, 0.0]),
array([16.13102990119284, 18.185565792692973, 17.54407926678902, 17.488980022335014, 16.734205730759808]),
array([16.291375336337854, 18.594840612327115, 14.82608633539149, 14.22488593249492, 14.800804063610238, 14.301305085617278, 15.5737810113726, 15.633783325224961, 14.04765413784575, 12.291159334933175, 12.853990095357595, 13.148190487264387]),
array([15.745447370903385, 15.211270153770375]),
array([16.247196513825436, 14.908553943300237]),
array([16.883110670338738]),
array([15.021755306108929, 6.275345848864566, 1.1745776004806365, 0.0]),
array([13.864759061017129, 13.797772018710177, 9.789597061399066, 7.410713123261705, 8.49064741021577, 11.511203407666116]),
array([14.82808779819055, 13.95497348527012, 12.601300552482176, 9.110544114756493, 7.7003423499984915, 5.847629687098566, 4.112412263869645, 3.790041848618243, 3.0355034959310663, 1.3375499243618667, 0.0]),
array([15.595946638213647, 14.68781434969287, 14.163075310472008, 13.404706748647001, 12.729895661973655, 12.746592645457714, 13.281059714220238, 9.992374780127442, 9.437855236486463, 11.431306773295264, 10.426369180887711, 10.073825407753029]),
array([15.331905084017404]),
array([14.859032063154913, 14.278760740627602, 15.326368492677998]),
array([14.694093486571138, 14.735401761518375, 14.096085631482481, 13.278870427477434, 13.158405482606055, 11.30445209220457]),
array([14.3953815575889, 13.892740260238533]),
array([14.370015072692484]),
array([14.700880900023297, 10.694704816224776, 8.469002062540893, 8.407334123398218, 8.298527605118325, 5.699809915114304, 3.7879072308445325, 4.161239827783132, 5.179078203896337, 0.8015988651509695, 1.7255250742006243, 0.0]),
array([14.262674905434269, 12.354608877814751, 11.661113279797998, 9.319192670897039, 8.792493916043892, 8.625509425072305, 10.638152168411425, 5.766304067665477, 5.671880036478527, 6.990057938455921, 2.8153318999871995, 2.204464325762535, 1.3589615615785666, 1.582670166360031, 0.913062742012954, 0.0]),
array([12.87585468414719, 9.269991008931814]),
array([12.24777665056477, 11.168116175198259, 10.599063673828033]),
array([11.946660581586944, 11.887431743426518, 8.526355844309025, 8.711321887240505, 10.697236682689079]),
array([10.432318705344208, 10.932481139755167, 10.370573127402233, 11.553376569041358, 6.76518051672122, 6.872667097864021, 6.822624737978075, 7.181883599899473, 7.179814172407531, 7.053654873780068]),
array([12.071091054709877, 11.06354398114702]),
array([12.740625522543732, 12.34822034025876, 12.570394136973217, 13.059653206678309, 11.775251280106332, 12.50617915464539, 11.395068517794646, 10.72171378903892]),
array([12.993905982162385]),
array([10.249301410063003]),
array([11.830691278262321, 9.951589362773435, 9.773618762918426]),
array([11.748015597773557, 10.574919308306963]),
array([7.5679333150077595, 1.4993818627456537, 0.4583775861978226, 0.0]),
array([12.138077083862985, 8.089432051654796, 8.697771833195672, 11.539469163239222, 9.718106428574188, 11.027003468956245, 9.534342073144783, 6.420624441884377, 5.91686262742781, 4.276395226182783, 4.706424607247562, 0.0]),
array([9.003631863776715, 9.43101157155819, 8.951252443208304]),
array([11.531305786683681, 11.457768437734712]),
array([11.768148221560796, 10.804521254390274]),
array([11.55062792862406, 11.253394294631835, 6.98474400577531]),
array([12.219272773683503]),
array([11.949304783439345, 10.141948984418065, 8.75637682411084, 10.731057015156372]),
array([7.692427678467078, 11.283134019014893, 10.345180595641088, 7.272279440303874, 8.74471889344081, 9.340990983111068, 8.169611660396205, 7.225771560821212, 6.620314119523835, 5.657226817330042, 7.133347909473124, 6.297908684969738, 6.955108320970295, 7.076079827823609, 6.708661362979255, 6.961933573331024, 6.787895243558062, 6.211650233332239, 6.879629270557317, 5.6400870760235176, 6.596879784479512, 6.575924206915989, 4.639269816458329, 4.954435469494433, 5.0467177618299734, 4.918462967748534, 4.986417579265205, 4.026881483777151, 3.505776060668552, 2.9587147638777593, 1.912777370421988, 1.6288169382388178, 1.579825142002699, 1.5276733855106057, 1.2352583296951871, 1.2629135029167737, 1.139386379792279, 0.8630944092512751, 0.8510212677675794, 1.3379078681810692, 1.5963218424877796, 1.3476878464701645, 0.6872188730170502, 0.4910655967666004, 0.07128427812260739, 0.0]),
array([8.684624040202298]),
array([11.592091745661081]),
array([7.739756547779208, 8.397877846579698, 11.130911846429642, 10.69751474737663, 7.3418591887275015, 6.22485112683076, 7.0345488418393485, 6.785985706186987, 4.343108885189874, 3.771767955619173, 4.634525843556016, 3.702800982214326, 3.858040774445138, 2.9854780823474343]),
array([10.563826696536221, 9.960095337907843, 10.279957804006976]),
array([11.18798361622708, 7.865602371592534, 11.121915975182679, 7.869454222591132, 9.310618738434925, 6.421808139718425, 6.32576500726042]),
array([10.222876346892217, 8.396686740170654, 8.872705725258037, 9.872271938743287, 10.45010679303444, 9.49554059554605, 8.208570994010028, 8.646629223259637, 5.55638202541331, 5.37997335660466, 5.862541801859185, 5.64563312937689, 7.029328484219823, 6.536500399230893, 6.769504701483559, 6.670338718801918, 4.5327220162830075, 4.172385921124362, 4.675657705395845]),
array([7.390284545106363, 9.539732522092942, 5.992632863210643, 6.828371476660405, 6.78167030877143, 6.405334814701163, 6.2179447792495655]),
array([10.402321052293187, 5.771356643165274, 6.801562874137151, 2.8738754276237857, 2.045304401459453]),
array([8.70882893398139, 8.262476206810755, 7.313562881362167, 7.534774247873659, 8.93849698826737, 8.294782005626839, 8.262021276904326, 6.203656950388716, 6.855325167479909, 6.402304690374327, 6.328189843861329, 6.645248415941871, 6.319309519954318, 4.909628603869605, 4.752626252307505, 4.683723334804377, 4.8895928458769236, 2.6679016758107212, 2.2809249841444186]),
array([7.346202280441649, 10.496047639158524, 10.464626347124902, 6.117829898577382, 5.675024794004748, 7.125380554326756, 6.367297216371635, 6.810390054975058, 7.0450675650892425, 5.4546288163608585, 6.411281993854107, 7.153435141189779, 6.703960702498421, 3.255139296162623, 2.0246592489778754, 1.755557733297003, 0.9547483461436855, 1.0382403342410158, 1.255473325611292, 1.4302886049823882, 0.5136441970013775, 0.0]),
array([7.4833465484698465, 6.764954617208367, 5.369701240140803, 3.6556185398059937, 3.0522109504320083, 2.8997231470139107, 2.2350422868992155, 1.5625363394573388, 0.0]),
array([10.693545851410454]),
array([9.696723870076115, 10.530672165257897]),
array([10.103993054835055]),
array([9.883184953721198, 8.51698040164922, 8.636114142404512, 9.493287515576581, 10.257168842052971, 9.98222445474336, 8.071478404593137, 8.147472275239863, 7.538257152213475, 7.38392307185396, 10.280699769889102, 6.454868935400067, 5.764015287325454, 6.764993758134337, 6.622356183767827, 6.306577355565383, 6.607804584891135, 5.909685332939808, 6.2367439380696155, 6.4820138809162655, 7.096885142362419, 5.728541772011967]),
array([6.2470357096580855, 6.282372251584842, 2.7660224027717284, 1.3539289499012086, 1.2404177679214854, 0.8556414795006605, 0.0]),
array([8.979315366673793, 8.950292562064753]),
array([9.565503708147824, 8.67374957363301, 7.262814888579497, 7.008312323801585, 6.37519619675761, 6.632717935264181, 6.816864216138495, 5.600925533134464, 7.0995570493190945, 6.645383220190122, 3.8205086427271016, 5.115714533168076, 5.262023347335649, 4.623561962502052, 3.5952668883211962, 2.6809403712492514, 3.4522751630235384, 3.555325478101379, 2.9849855213844694, 0.97798916943065, 1.783190647208438, 1.1298661362737046, 1.6453330931718924, 1.4288761014632325, 1.2655098008911756, 1.6104642997695773, 1.0981461135693644, 1.2591728182096238]),
array([7.993955949839674, 9.165717330719277, 9.152716350416846, 9.411041217777482, 7.69361444454505, 8.184635616603607, 9.051372401270937, 8.137138910264825, 7.562156685479721, 8.526435444839574]),
array([8.86566577025753]),
array([8.575990457985375, 8.380120787544]),
array([9.131423557864638]),
array([9.025782715742318, 9.11236294760348, 5.053358276851904]),
array([9.038168331141716, 8.711529267056845, 7.122918142961757, 6.607315104124199]),
array([7.500731052109026, 8.693268420519514, 7.041345703430451, 6.940747195303643, 4.635623743015168, 3.376567133250484, 3.3602009578238246, 0.8047676826625209, 1.035758223064661, 1.2369463248207548, 1.7502886169858514, 1.4663483734535472, 1.1307685660508071, 0.0]),
array([6.126677290639418, 3.1197564238906095, 3.205613945925411, 3.3920596831954715, 1.5818679117782732, 1.284515776087214, 0.8486971960329327, 0.0]),
array([8.017138812575054]),
array([7.3779717385001105, 7.988902577008702, 8.050993452695366, 8.542788639364803, 5.6139805227089665, 5.648885977055114, 7.033071402509179, 6.262389223457804, 6.987869794542249, 6.090469591607634, 5.569718902924599, 6.971914845054524, 6.554560254813571, 4.509391999012397, 5.2792532650533195, 4.4374915897365685, 3.2943527013040237, 3.130505957622866, 3.4318610870748674]),
array([8.593942858480498, 6.877747301105963, 5.585613318858691, 7.02348520376238, 5.759632440198387, 6.872842798059852, 6.617434918983892, 5.974081626129739, 4.208373099363454]),
array([8.675737516825974, 8.42387813333927, 8.858636666022242, 6.2843820717015, 7.152603679488133, 6.3567753063925, 5.651755604972086, 5.9724350346144455, 5.393357597030161, 5.884940644458892, 4.831947523796943, 1.2271539511298502, 0.0]),
array([8.946170713218047]),
array([7.61122301066797, 5.80051591652853, 6.463563072613322]),
array([6.815753619198786, 6.538067448148479]),
array([8.548105035944747, 7.53434350445441, 6.375459680761377, 5.762109954219667, 5.809450661917085, 6.582212071588974, 5.471110230129334, 6.744615958426691, 6.254762186895908]),
array([3.4750676351380543]),
array([8.014723768421849]),
array([8.56369912262881, 6.57212021783624, 7.20716696872162, 6.480951905302983]),
array([5.806109085581926, 6.633962919753334, 7.12430911892521, 5.410479905666587, 5.535450937557812, 4.823797061276549, 3.7869162191489623, 5.196493849334898, 3.085397898475027, 2.123085261227931, 2.024434310754251, 1.4317587240640433, 1.7249792782581657, 1.323257974492193, 0.5278414855805819, 0.24400253237768976]),
array([7.56222008987535, 7.255353222004601, 7.6240499477714705, 6.506984367288669, 6.3087007027730335, 6.420217808591872, 6.953089364409532]),
array([7.680669532735987, 7.2749753905942525, 7.8835739948811225, 7.8567745466436385, 6.13929343604603, 6.817288023156418, 5.511499238875677, 6.554258152972834, 6.513491354108356, 7.2384179276377445, 6.754133688355796, 6.372199102930643, 6.4478046710845325, 6.68992497484661, 6.870933805908068, 6.678120381889398, 5.825086344152313, 7.070206849513946, 6.245064613247812, 6.874646255560407, 6.702517164764098, 5.6266379884776025, 5.984413611658235, 5.808067808339143, 7.20941084138252]),
array([5.61645362546132, 5.984825756939283, 5.846653995663449, 6.902523496407507, 6.697401684541641, 6.1820653360103535, 5.634834848303557, 6.308158429313823, 4.42692477990679, 3.062303108739232, 2.945233209085356, 3.1895774528582774, 1.2980621990455048, 1.7643939383152474, 0.8327371170123459, 1.2696278028125327, 1.4838808665349534, 0.8807446473313201, 1.0433172143522493, 1.2408940475449233, 1.5945533568162178, 1.3237790237550082, 0.6305621139804412, 0.140681953556273, 0.5653141820147933, 0.0]),
array([6.345897500071457, 5.9520031999812, 6.0573879122021195, 6.71462316619651, 5.33910584497393, 6.231788038544293, 5.775958171654077, 4.118482476615437, 5.237207958631799, 5.089646132170908, 4.0581861029442, 4.82208440647316, 5.296793078487501, 4.590201098402055, 4.551725717563805, 2.8803354549013287]),
array([6.439654594369859, 6.841121200439103, 5.451956141906093, 5.876643377005584]),
array([6.238920590110808, 6.58249521931123]),
array([6.878565453866292]),
array([7.316253545518296, 6.8759255013599025, 7.063280211269893, 5.786947010092847, 5.583250324095146, 7.097338989134988, 6.45637905727244, 6.135855311365078, 5.862746162087779, 4.48967959705868, 5.102712420547888, 2.6375524529324688, 2.944262506336684, 3.1219651708170577, 2.6131992911210156, 2.6614826938088334, 0.9341579646746694, 1.7321146118645852, 1.359482120584494, 1.184843097285908, 1.4163703028583905, 1.4045463021864952, 1.072138075518264, 0.5123098189371762, 0.025939953024060103, 0.0]),
array([6.524916686426837, 6.35384848706275]),
array([6.884645644344217]),
array([6.225093212422479, 6.862794872178766, 6.750856789809392]),
array([5.514291448123259, 7.0777741046146545, 6.367643015841436, 5.416524592405906, 5.490661017725547, 4.779363245284253, 4.202111463587201]),
array([6.941733240565833, 5.5831033871562115, 6.142444106073283, 7.206991857097448, 5.359266298810255, 6.996392997594281, 7.163529885914973, 5.574533643475512, 7.235754775719515, 6.843961517441241, 7.026929647728708, 5.542720289053232, 5.680105345675541, 7.158982385052194, 6.99569498277112, 5.55690485929202, 6.30964858013627, 6.538699703162289, 5.192427190302641, 5.268405392986808]),
array([6.805677074236429, 5.58493179377067, 6.251376131406364, 6.234684965591327, 6.462200345615177, 4.137451883758436, 2.671995934273785, 3.4408081824866072, 3.0290839834429297, 1.5746673152647452]),
array([5.379283798133069, 6.600713309476173, 5.7356428366066705]),
array([6.71704993034279, 5.733153471653765, 6.633970519016797, 6.641427695578855, 4.907565104638379, 3.4427689935555197, 3.5125572964775715]),
array([6.321640994324716, 5.596974628277815, 5.802082738115562, 6.7612628994923085, 5.737300509097714, 4.801181814069161, 4.14426215990026, 4.3594736075764535, 4.491476946006238, 5.005325081961802, 3.8618047012294068, 3.57081136163052, 2.8844741766708863, 2.9386672534861757, 3.481744356748264, 2.3355602034623266, 1.8127313648143482, 1.8490920248907496, 1.0928942020833723, 0.8932155643450652, 1.1396865190594005, 1.7761167547852736, 1.184369401678138, 0.23400766693663, 0.0]),
array([6.297127847470931, 6.35866968788686, 6.599927380087164, 3.7883806840893284, 4.77439534509709, 4.294536756924609, 3.8559729490811794, 3.1656658545176803, 1.3084156193827097, 1.2777096844856786, 0.9021524023417281, 1.7608266463046525, 0.0]),
array([6.627880343654019]),
array([5.376246727633719, 6.519115339180058, 5.573186575564307, 5.333272926472707, 4.680185797915561, 3.834156881090621, 3.22640942377373, 3.3974461003735783, 2.8821225965285655, 1.9300118152041357, 1.7585854511387258]),
array([6.224642941326756, 5.369172695617208, 5.303297561715367]),
array([3.505026261117702, 2.867130370381202, 1.6605786469075303, 1.561886554059212, 0.8480322733004128, 1.6135264254267974]),
array([4.556203110932956]),
array([6.136480821509418, 5.969133811824226, 5.351621714924033, 6.228515982158668, 5.34822028506482, 6.359750144532948, 6.253460831545242, 4.083795683684531, 3.721803754059448, 4.604184053888657, 3.7561290743103406, 3.442858794991028, 3.282147227305049, 3.4256906717274966]),
array([3.1178025796654025]),
array([5.936698004233261, 5.430494961545863, 5.433282429944756, 6.189143156678795, 4.2717379626034235, 2.6397212947652426, 2.9550010928864747, 3.0923814642293106, 2.4436989809230094, 1.912826991433223, 1.868303770556251, 2.57780762188348, 1.4921571145018293, 0.9522940964301398, 0.8597839785176724, 1.4286337218938439, 0.13469392421406634, 0.5959257383375977, 0.0]),
array([5.944635980279442]),
array([1.3236480488056985, 0.0]),
array([5.4426254498010875, 5.438362543087265, 4.720852884699815, 4.707308346617783, 4.395643490538151, 4.036654802370528, 4.526589635262273, 4.047483538144494, 3.0838014511956584, 3.118738860404603, 2.9639180879569573, 3.1217482492954667, 1.1852298928061131, 1.2422298798610096, 1.6176439411235541, 1.0138622766507637, 1.536438415759538, 1.7859676417149188, 0.9917735976418459, 0.6823989564047136, 0.12868773927584487, 0.5610308449201238, 0.22482421750426307, 0.0]),
array([5.807743937223821, 4.069581978527619, 2.6627310654775815, 2.681915076879487, 3.06213716090957, 2.7030289065813324, 0.9551003156025104, 1.1999858431755597, 1.6155506303797444, 1.480129498698643, 1.6759264437750379, 1.7321425746117833, 0.0]),
array([5.56352285429481]),
array([5.600230004950311, 5.340160734887382, 5.4261770315400035, 3.970173510079973, 4.1391995979438825, 4.573780790838869, 5.09361942574241, 5.236432020097843, 3.477865146210509, 2.6450508682536507, 3.571624170842424, 3.0917417705414003, 1.788083949903383, 0.782493941440668, 0.8896197894669196, 0.9953631020541877, 1.7813671815005576, 1.073086678962622, 0.9449932945168897, 0.0]),
array([4.876087164017671, 5.120267563188944, 4.378922766789483]),
array([4.394718291841723, 4.689901555124392, 5.175006593448025, 4.53645393819207, 4.587414955371323, 3.149068962475056, 2.6297009628074677, 3.3270002446040965, 2.5976129018910314, 2.578467577131564, 1.493460913349681, 1.6203966403853485]),
array([5.627814479950421, 5.655472231412002, 5.328568671944372, 4.947686229505452, 2.606819777215937, 3.0075441164027246, 3.405520383482363]),
array([5.732569105777786, 2.824806609211416, 2.999215961294946, 2.7087266351663692, 1.173825758471585, 1.3711454044554574, 0.8256333860489311, 0.9007217983177611, 1.0099453891242653, 1.6337487777211854, 0.0]),
array([2.836658742723662, 2.843827471718773]),
array([4.969739798163515, 5.217497893104583]),
array([4.21552498970634, 3.341248926759071, 3.5114881219455074, 3.172910037890893, 2.251374622017296, 1.5925782231373398, 1.3066648257048574, 0.0]),
array([3.9714054784579718, 4.235212108351007, 5.1686032373381385, 5.187477790178871, 2.5948109983931094, 3.4905524566494375, 3.1419955390295438, 3.5383060494509535, 2.5106872556567965, 0.9588360747853969, 1.4363652752173197, 1.3548042350249911, 1.426481975807103, 1.351917142487146, 1.3323500884534227, 0.7083082304796363, 0.654660327637144, 0.4586973451963703, 0.0]),
array([4.277846116565058, 2.8083142261989904]),
array([0.8442710749527474, 0.0]),
array([4.25179689612292, 3.9178310910716805, 4.035730228465431, 4.607423669309991, 4.826757085000173]),
array([2.866233419464935, 2.7717946084008958, 1.131522595022989, 1.2704907478449048]),
array([4.041070037517647, 4.575329955669677, 3.135381077251927, 3.517845633647544]),
array([3.674709826287085, 1.3332472084912386, 0.9184927546305516]),
array([3.447730698008156, 3.5003150948637463, 1.2336258265779474, 0.9061232935723269, 1.6641400236338302, 1.4113272208038345]),
array([4.25570277934708, 4.462402380961927, 2.974308848341539, 2.6909658533666434, 2.141110645636526, 1.807308808517205, 1.4384914770279584, 1.7356931647758485, 1.7758060266992615, 1.187374347053091, 1.2824289981694978, 0.9508500846286143, 0.8145842404260509, 1.5203703715330068, 1.668881809050244, 1.548899743786188, 1.106762480395261, 0.40484170918558104]),
array([2.9087330535715124, 1.0543921937367018, 1.7448922997390675, 1.652520703489797, 1.6936809659378014, 0.23655223951974902, 0.0]),
array([4.417179341974848, 3.896726992828774, 3.7109279798805037, 4.397850009273051, 3.6431407503964204, 3.097596759685799, 2.9132064139512637, 3.209821564245205, 2.6145534830514605, 3.51009111024224, 3.130300229979834, 1.2510731249678613, 1.2032713467137732, 1.7090239954118525, 0.8179842087528599, 0.9443045712603975, 1.6283816387317767, 0.9473179019809106, 1.0281889837962461]),
array([4.400990580609069]),
array([4.396904352660192, 2.767481503047524, 3.252585617074813, 3.383444456406451, 3.1376892723420227, 3.512553651182966, 2.1157144519520426]),
array([4.201586843697646, 3.8059423183366117]),
array([3.792101354599989, 4.05315049094468, 3.49396394391523, 0.9253931856897978, 1.6780270720376995, 1.4654405464952642, 1.0251048324987966, 1.7898462864450027, 0.961601847518133, 1.4800604092754193, 0.47776747486075005, 0.4807584651781419, 0.0]),
array([3.931984154912065, 2.924689226812645, 2.719588795242788, 3.5372075924610593, 2.3831739067950086, 1.4460343584791135, 1.2811416091615575, 1.5636729044180453, 1.5980595019348098, 1.7886703183677781, 1.1557842978708972, 1.3084087555039243, 1.1891328634182916, 1.4313122759354635, 1.3690792176437796]),
array([4.10301377560022]),
array([4.018613951656964, 1.0837109994304857, 1.4899165278898998, 0.0]),
array([3.1988840300385197, 3.062891845482459, 2.6797991758762967, 2.0457456487819643, 0.9931493210731129, 0.8032384869693858, 1.2023987655775958, 1.5477573591393594, 1.6715450899176618, 1.080943074769845, 1.6547326456396712, 1.0022576897292041, 1.0336192588441784, 0.0]),
array([2.962327603176241, 2.393998621808997, 1.7277290965681371, 0.8420598283563873, 1.65206499179142, 0.0]),
array([3.0078493493249914, 3.194987695727616, 2.9055446647805283, 3.2959995943149796, 2.8167785894168946, 2.956994749789303, 2.953533342418236, 3.057483021097732]),
array([3.704378554884124, 2.832565576086857, 1.0700771142568741, 1.313746006574521, 1.0240060778052962, 0.0]),
array([3.523992296901045, 2.715547526918185, 2.6609969548930184, 1.7293767936921711, 1.0527304870325205, 0.8110606869062242, 1.6316981101838457, 1.681664885457974, 0.9728345176312343, 0.5767320049476202, 0.0]),
array([3.1807067396517463, 1.1298517012218563, 0.0]),
array([3.46807620241542, 3.239589540329383, 1.66386393300355]),
array([2.4898275018169755, 0.0]),
array([3.336379457171605]),
array([3.021133072797788, 3.107145634488675, 2.9136378324309953]),
array([2.4021153034250355]),
array([2.8506268066183025, 2.192846004300053, 1.7658620165616485, 0.5496097967653375, 0.0]),
array([0.898497607417563, 1.4950239736446491, 0.982562185114208, 0.4136250319932611, 0.0]),
array([0.8354756305608853, 0.8118904782585648, 1.5447638291731711, 1.7304008790522787, 0.0]),
array([1.4841561388588667, 0.78346833523363, 0.9675029088157606, 0.9661878400603385, 1.5019575580735232, 0.7851378429004774, 1.143200518074953, 1.5766469454024847, 0.6569487466611026, 0.0]),
array([1.81743260470402, 2.28466769198759, 0.0]),
array([1.144872390255415, 0.0]),
array([2.454642765909461, 1.5888500832095886, 1.616752502927049, 1.0321793989644177]),
array([2.347892900349151, 2.4722684506200965, 1.321006667431475, 1.4093945821639424, 0.886527327581455, 1.6160991439525905, 1.5389497375994694, 1.4914292882740068, 1.1513661262594117, 1.2224958818300566, 1.395422948435978, 1.2236601869153243, 1.0456606898251706, 1.499843762395246, 1.1894687924107616, 0.7780621194446042, 0.49660499126993074, 0.6369791072248546, 0.0]),
array([1.297427911024565, 1.344655226359496, 0.3677006050603042, 0.0]),
array([1.7274751628946359, 0.8547671464508809, 0.9687920902383884]),
array([1.8185216316881652, 2.0426456470794343, 0.9751583114306295, 1.4796503851390939, 0.0]),
array([2.3594674735986474, 0.8997411686743872, 0.3828612447405552]),
array([1.2037800445219964, 1.4831241483233943, 1.433226139321039, 1.3705764619273473]),
array([1.0406892928799474, 1.0219381539693542, 1.269955922855976, 1.612089608674001, 1.7353246402232, 1.0192034694070626, 0.0]),
array([1.7858904426748683, 1.7182423627696986, 1.435040767412857, 1.4278371532162712, 1.1872290490804822, 1.4885002598307562, 0.505765816786887, 0.06010694119318914, 0.0]),
array([0.9598367147951653, 1.3917775648800736, 0.0]),
array([1.644248039636476, 1.1019143497228998, 1.2154594861514343, 0.0]),
array([1.9365361221044277, 1.5590745620001032, 1.359137151069393, 1.125798457767826, 1.2243525166133002, 1.3473621650984564, 1.0674996378531496, 1.1757021233928202, 1.3664077684899112, 1.263689037586891, 1.5859953736878054, 1.2752467945890031, 1.4127227171956447, 1.3128587405504581, 1.2378738631150827, 1.6973435398868808, 1.1947554382347612, 0.24332089358129594, 0.4845510057469962, 0.6256354310883141, 0.3777078727063697, 0.6123577304284754, 0.5790547704850119, 0.7497140356710229, 0.0]),
array([1.5052180946342615, 1.6978183519653374, 1.3189247030215192, 1.6974335342631055, 0.0]),
array([0.025291945775266098, 0.0]),
array([1.4244589409500545, 1.3986477151234609, 1.6118903290601319, 1.5650762665524862, 0.8227646619763143, 1.4199875234042922, 1.7518991212684427, 1.1375460553379901, 0.900443729695191, 0.5484588618882964, 0.4778940741335914, 0.0]),
array([0.9954793746018006, 0.9745142793974797, 0.0]),
array([0.78300637623072, 1.412824164280179, 0.8530925641938832, 0.0]),
array([1.0043805298377726, 1.2474808926861103, 1.1111819205769273, 1.0281208306279932, 0.9574331731985235, 1.6961525580571373, 0.8809277020946351, 1.0771339076741226, 1.6872161615372272, 1.07883378578655, 1.4234032338708464, 1.4250944987441188, 0.5567708503478355, 0.10133196788440002, 0.0]),
array([1.5917984612736424, 0.0]),
array([1.5052632615026134, 0.0]),
array([1.3876617216783633, 0.0]),
array([1.5168228865376796, 0.7908889349394895, 0.8221805616005158, 1.0417702437228316, 0.7982353442020855, 0.0]),
array([1.4235876183691, 1.3822894015998153]),
array([1.3975414372348927]),
array([1.1959816834156334, 1.1538377641000472]),
array([1.3557982456587405, 1.1699601941469653, 0.0]),
array([1.4087079553173956, 0.0]),
array([1.1141887587161097, 0.9548805682468213, 0.0]),
array([0.6755790666863264, 0.5538890648874001, 0.0]),
array([1.1352602552975029, 1.3416724091835646, 1.2521520017934875, 0.0]),
array([0.9685160882067574, 0.9780819546946844, 0.0]),
array([0.8608731813380872, 0.7292325907649089, 0.06949288111652223, 0.0]),
array([1.0442031973955856, 0.0]),
array([0.8321269109144334, 0.36391993691362756, 0.6995851698734672, 0.0]),
array([1.047117169316367, 1.2203137304502625, 0.0]),
array([1.0510045308967637]),
array([0.9879108700424383, 0.0]),
array([0.8018432355733291, 0.9492625198948371, 0.0]),
array([0.8143962324843399, 0.794894320186895, 0.07986941574267784, 0.0]),
array([0.8034486574664075, 0.7587878197356281, 0.0]),
array([0.5748945707250626]),
array([0.6442680408478887, 0.4763656481084669, 0.0])
]
d = [data_1]
names = ["21"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T1', 'T2', 'T5', 'T7', 'T8', 'T10', 'T12', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T25', 'T26', 'T27', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T40', 'T41', 'T43', 'T44', 'T45', 'T46', 'T47', 'T49', 'T51', 'T54', 'T55', 'T56', 'T57', 'T58', 'T59', 'T60', 'T62', 'T63', 'T64', 'T65', 'T67', 'T69', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T77', 'T78', 'T79', 'T80', 'T81', 'T83', 'T85', 'T86', 'T88', 'T90', 'T92', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T100', 'T101', 'T102', 'T103', 'T104', 'T105', 'T107', 'T108', 'T109', 'T111', 'T112', 'T114', 'T115', 'T116', 'T117', 'T118', 'T119', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T129', 'T130', 'T131', 'T132', 'T133', 'T134', 'T135', 'T136', 'T138', 'T140', 'T141', 'T142', 'T143', 'T144', 'T146', 'T147', 'T149', 'T150', 'T151', 'T152', 'T153', 'T155', 'T156', 'T157', 'T158', 'T159', 'T160', 'T161', 'T162', 'T165', 'T167', 'T168', 'T169', 'T170', 'T171', 'T172', 'T173', 'T175', 'T176', 'T177', 'T179', 'T180', 'T181', 'T183', 'T186', 'T189', 'T191', 'T194', 'T195', 'T199', 'T201', 'T202', 'T203', 'T206', 'T207', 'T210', 'T211', 'T212', 'T213', 'T214', 'T216', 'T217', 'T218', 'T220', 'T221', 'T223', 'T224', 'T225', 'T226', 'T227', 'T230', 'T231', 'T232', 'T234', 'T235', 'T237', 'T238', 'T239', 'T240', 'T241', 'T242', 'T243', 'T244', 'T245', 'T246', 'T247', 'T249', 'T250', 'T251', 'T252', 'T254', 'T258', 'T259', 'T260', 'T261', 'T262', 'T265', 'T267', 'T268', 'T269', 'T270', 'T271', 'T272']
def get_taxa_names(): return taxa_names