#!/usr/bin/env python
from numpy import *
data_1 = [
array([33.5943410757728, 33.053475778809386, 33.454294741935755, 34.3152185301979, 32.789092394025474, 32.83936169656261, 33.09016461918195, 34.511319143899065, 34.59478390002796, 34.049965035038085]),
array([32.44208477657085, 34.42850892396132, 32.0808260486778, 33.45236540062682, 31.20302969304816, 33.2611979732592, 33.89810776883316, 33.15289314810913, 31.131005533146592, 33.22485704141337, 32.36360320575659, 33.408891123020645, 32.56157469485884, 34.55714999753463, 34.52250259250486, 32.86323161125329, 32.38853631845555, 32.30803838611074, 34.24099547250638, 32.08476061775994, 31.975892585685077, 32.98448115893937, 33.28609624595253, 31.47711536226848]),
array([32.48197571760377, 31.930388688406378, 32.81588366784957, 33.47406981815049, 33.28489308352972, 33.081187077674265, 32.73955477113157]),
array([31.046875936870787, 31.85836163975387, 30.267618909608174, 32.94268788099388, 31.32555267807847, 33.22222175229435, 30.63654727407977, 32.48202255303111]),
array([33.81044192319497, 30.760189789704445, 30.857137621351804, 30.32353969733787, 33.16487274318226, 29.79227010924896, 30.719069283669295, 29.271481266156798, 31.569741374104282, 33.302190135082256, 30.37421952012743, 28.529373174293113, 32.92631058163856, 31.922103176898382, 28.288535554790887, 31.488471389955237, 30.315019938754684, 29.82382363010729, 29.267401315453426, 30.169155868060063, 31.747973089918503, 31.89741790229157, 29.88128251839204, 31.039647325482786, 29.676092722309775, 32.223685832007654, 33.349723804989324, 30.10293783014885, 28.103724438003955, 31.057473338559227, 31.43927064400269, 32.6517459142909, 33.30549577627589, 32.5456041308849, 26.47081167906337, 23.9044523973216, 27.41365038686573, 25.543240689299047, 26.19370220363775, 24.53338644102836, 25.01722191977779, 25.818018655996184, 26.793536402867545, 24.06193554408224, 25.280668914053972, 23.591592798089437, 26.923797466748614, 25.665058301801054, 24.72981947069842, 27.566573880424496, 24.754563913988388, 22.19009121392755, 21.34124215566535, 21.458901855650073, 20.974858211036306, 22.521081844584458, 21.921053301632714, 22.820546066875707, 22.744296373396892, 22.141558938576804, 22.59558340784827, 21.644725811488627, 21.433947973779667]),
array([29.94243971195912, 30.449889097844885, 33.21765219294947, 29.959484456414923, 32.16198123686878, 32.07110891769729, 29.03767675210806, 29.36022698230838, 32.74559973632741]),
array([31.283874705063955, 31.66630473334586, 32.488408119934576]),
array([32.016600996874, 31.909274061280218, 32.17109583700035, 32.0644698858828]),
array([29.082306969195812, 30.597447921936595]),
array([29.25620540127958, 29.49857230236808, 28.10482140247745, 29.597620699926694, 29.63122605814586, 29.160924511235844, 29.62368821645872, 28.575515663356214, 28.217001102086755, 29.663984961636952, 28.33228215540151, 29.167912835850956, 29.006672373206943, 28.588226447956853, 28.341077484630738, 29.343808001819134, 28.84498289076345, 29.09083737972866, 28.764801434044877, 29.585922017656735, 29.43504767404662, 28.190320299207734, 27.93092635513467, 27.427205150727698, 27.626801908703268]),
array([28.62202734274424, 29.38010234945713, 26.56174466744272, 27.701810836241403, 27.566778035370195]),
array([28.557415129679143, 28.22641440153901, 28.682927269638085, 28.46890643638101, 28.541690670999277, 28.134094463162, 28.373317475114636, 28.16324154463272, 28.614641857156148, 28.724597214546947, 27.83615480431569, 26.26912512183915, 26.36190822586209, 26.200531611381006, 25.57010804445401, 27.503922337483267, 26.637933623886642, 25.753091588740947, 27.582426289703136, 27.35462214410569]),
array([28.735604725964347, 28.691670696205442]),
array([23.766021641650003, 25.44558289849129, 26.428202046812192, 25.74745267616084, 26.259964438523372, 23.26096535326602, 25.37192875225893, 22.69017522763068]),
array([28.661733889249714]),
array([28.02715754969639, 28.055500769674023]),
array([24.494404313216243, 24.495061095957205]),
array([26.199619222654178, 26.71818169104887, 26.314365465600176, 26.622527920135465, 26.376897439971156]),
array([23.681737613356542, 23.923405005841303, 24.443514384750117, 21.480607821527308, 22.618924294673068, 21.966155825879973, 22.065249546452, 16.103883776629193]),
array([24.285792387759845, 24.101944215596237, 22.255455079371934, 20.718979657064224, 20.614907196921138, 15.12217241156187, 13.195875791254714]),
array([23.796465626721044, 23.242397220290407, 23.895238368623904, 24.29532777295084, 24.11223008576703, 21.964516898169034, 22.99438493179766, 22.853043067250844, 22.585429475335655, 22.21577964235153]),
array([23.471464466684534]),
array([23.600887532573225, 25.3042348214499, 23.633817947477873, 23.92615059603859, 23.971444678099683, 25.25547748542358, 24.23325658189855, 25.22351997276194, 23.821422396568373, 23.502022399228075, 25.13719169501848, 24.43866082686111, 23.839106533649304, 24.76014610482833]),
array([24.802380830070593, 23.90929907973328, 23.49918850503066, 23.730713761015362]),
array([22.31880290818341, 21.91192372922839, 22.56691302431507, 22.0917689599815, 16.721787662145662, 18.27242598417127, 14.809868162265085, 12.589698476756015, 12.51746008537701]),
array([23.479646679155206]),
array([23.035678027745046, 23.093057497043358, 21.712162685841193, 21.259661234800188, 22.361040743390205, 20.708296762782208, 22.818424464195743, 21.893571457889397, 21.995511552287578, 21.730245158785557, 22.573886317502296, 20.562785992093644, 22.737099254481475, 21.474558579027292, 22.726630399807995, 21.31461291056579, 22.386734068037953, 20.710184937958292, 22.923277341140647, 21.730513709057906, 20.500151072012773, 20.64560552991295, 22.584450416210505, 22.34015390524036, 21.77587662657841, 21.7769259331219, 21.084327577639463, 21.015190394372876, 22.452131296557535, 21.03196747245582, 21.218091596046552, 21.898164730110302, 21.453165375162694, 22.25640708381289, 18.197086182323247, 18.747939052810292, 20.348294173145526, 16.01114505125285, 19.74909224477392, 16.514273621477678, 17.749288035649418, 16.905579389743348, 15.20364801678952, 15.452204212367496, 12.914865680108338, 11.862707304841038, 13.228624728401789, 13.612336460968088, 12.32146502224353, 12.721726936625174, 13.267614597268048, 12.107217210627384, 12.130494439911368, 12.0755020620058, 12.783019160646319, 11.738571245127392, 9.168162213375354, 9.802774775694179, 10.892647174890223, 8.784801634454482, 9.326393398640306, 11.096537475534905]),
array([21.755962741149876, 21.753669812613005, 21.995613801877422, 22.224200725005566, 21.958970136883263, 22.34627427481245, 22.667263746836355, 22.775060123566995, 21.68271107172948, 22.159264061825784, 22.748233173628567, 22.08050831948724]),
array([22.75980181069185]),
array([21.3869534034253]),
array([20.874612918438224, 20.907866548964734]),
array([20.697405008475467, 21.676542362396106, 19.837258938799142]),
array([20.741755207504266, 20.473194413325146, 20.822490084206727, 20.961192774068607, 20.547676118560503, 20.649988147658917, 20.95918315041781, 20.90759352061584, 19.808450233385464]),
array([20.459999033433885, 20.590896307128894, 14.33727806362013, 12.297232447569801, 12.410810349691474, 13.461189463046221, 11.926193487172732, 13.692148303001261, 12.962513171016536, 8.34914274796268, 9.815911398222793, 10.471946130533164, 10.778255028134623, 7.495755313979269]),
array([20.455403164053923]),
array([20.157759289922033, 19.07343104043313, 20.248672974148196, 18.890796617252718]),
array([19.795920765363046, 16.16473280945832, 17.790185419855383, 18.822280864631942]),
array([16.60206864497022, 19.40570245900019, 14.329077884444464, 15.080373329792543, 15.02083330643397, 13.007045876047387, 12.740438227812254, 13.62114064721183, 13.60680627747229, 11.961188380596786, 13.405727919246308, 13.495772327843767, 13.611157497723806, 13.041451341531758, 8.972002860185952, 8.70379521202437, 8.23007215585608, 11.08518294047544, 9.049382752136307, 9.86992007918034, 8.398857052489731, 10.424542895535325]),
array([16.756235486962694]),
array([18.15089030908029]),
array([16.827337529837866, 16.348123796950492]),
array([17.65172821639631, 18.296794021232092]),
array([18.140448246779147]),
array([11.832576752211725, 0.0]),
array([13.182582597311558, 13.494976581697076, 12.955844197247076]),
array([16.664721387002896]),
array([15.997847458419848]),
array([16.17768159488004]),
array([14.843488536583617, 12.940058055150233, 12.571402034369106, 12.929074070107976]),
array([15.548853162446864, 15.247819438452192]),
array([14.565769927580225, 14.795545270563291, 13.713997167694256, 12.759121682833197, 12.78848930462657, 13.137129170905604, 13.32854763764403, 9.38841333943333, 9.887686849992752]),
array([11.763042830718339, 12.74849107443786, 12.133088797235448]),
array([15.495206113161728, 15.343297020776443, 12.236995214121954, 12.857650480195257, 12.101233519466449, 13.402529116998291, 11.908360664211994, 12.439589504901559, 12.483415207912495, 12.499945060557353, 12.80344312103217, 12.855711202286141, 11.867476797913257, 11.733626964460182, 12.342863789855283, 8.272475603361023, 7.834230507657653, 10.766617497385292, 8.6176907520786, 8.554273009958916, 10.382227634005279, 11.019238839153685, 10.845624055484823, 8.727980887049249]),
array([15.124531612670554]),
array([15.065437269829424, 15.20480702080363, 15.042641595927101]),
array([13.39646709831663, 13.42399297551609, 7.201831611507034, 0.0]),
array([14.13457745742599, 12.662940208870406, 12.76953680586329, 13.489284379881573, 12.427112482933135, 12.453117994230784, 11.353334546454105]),
array([13.039197193827013]),
array([14.053320688636928, 12.529428522041218, 12.979001979559191, 12.495726915984767, 12.763531431458839, 12.33522584181138, 11.793409038900483, 12.733627601619103, 12.51113593344167, 11.310030204952133, 10.968037024695075, 10.135537517038571]),
array([13.73777052189693, 13.686833049788158, 13.788446950269336]),
array([13.99360581072446, 13.388763008565627, 12.114663744811313, 12.312184001883649, 12.357006085091934, 13.334350145079968, 12.21383453657615, 11.83275223543395, 13.046851902938332, 11.856909936749547, 12.807961415260811, 12.653798157354418, 12.916491337790916, 13.470607521771784, 12.352622282588397, 12.056131576692215, 10.163187100555806, 8.088958946396644, 7.874514236385881, 10.85600183785166, 10.307601645301324, 9.836884155638593, 9.592639221214263, 7.8164435592515495, 8.58590813362, 7.22046411884545, 5.863371987061472, 6.672709678517144, 6.128186385815068, 6.802817415241222, 5.91549766181866]),
array([13.379554434811915, 12.650270574199245, 13.450859583574237, 12.092087656267642, 12.136028491197978, 12.979784743550315, 12.271453627373639, 12.390101846057373, 11.696215504267956, 13.718899982610022, 13.65641149643209, 13.756070152026918, 12.707907380725409, 12.783097169801593, 12.496537116669002, 13.38309474570285, 12.963108726717271, 11.786507528111343, 13.34220672578851, 11.34982300791274, 11.466201931716316]),
array([13.708431279465328, 11.948671795647748]),
array([13.17281312502124, 12.686916518575117, 13.400231073317062, 12.196798381681264, 12.877502467425343, 12.003350095816643, 13.28190450348047, 8.208649266781102, 10.662422652550864, 8.396143643397533, 10.70124070799164, 8.471606723099192]),
array([13.059555408960371, 13.267518255637743, 12.799709980872215, 13.5470439019485, 11.665019236962314, 13.549090898647313, 12.296480076704174, 13.028574728611986, 13.500867012746633, 9.113098986216679, 8.789047233645444, 7.3029741335623894, 9.28734984033019, 8.76837094406889, 10.494613165453421, 7.303623621388686, 10.226396603903625, 10.123021215014454, 7.374981838986503, 6.121317666600967, 6.318809768844187, 6.69016009676471, 4.250840107312021, 4.852995460324488, 3.665511140378426, 3.2068188060121807, 2.6536576513198016, 3.3870171022945934, 3.3724632088951214, 3.0334729275677295, 2.196603255340145, 2.437632566260548, 2.0470774954908983, 1.4044422392773552, 1.2471198870019182, 1.4611175119208302, 1.27657068570586, 1.5791300211915749, 1.4130823286257412]),
array([13.223623339080156, 13.333519454242621]),
array([11.801384444559032, 9.343915898224125, 6.656649393689232]),
array([8.080270007634752]),
array([11.460880039579203, 7.810243112819975, 5.571872408742472]),
array([12.195830038502093, 11.686963308063945, 8.509407907749202]),
array([11.726559455909095, 11.1229695459783, 10.873021810198424, 10.511089175895531]),
array([11.675841637130043, 11.667364212909089, 10.432232394650537]),
array([11.295909929082502, 10.049711368350712, 9.827618109357442, 8.509122618545724, 7.867853847033462, 11.284824576060993, 7.401958646062251, 9.437923782361821, 5.874605441368098, 6.35796138216805, 5.970162068812002, 6.974994629598537, 5.0466977220020075, 4.243289123750401, 4.974670509668581, 3.868060263813862]),
array([9.689305680892556, 9.273595460986039, 10.45674281705985, 7.3493257328540444, 10.062244154701629, 9.92956748512141, 8.996332642155467, 6.409315847853086, 5.597455549461551, 6.233131219876922, 6.651640346412193, 6.639536933139782, 5.323961462320997, 4.2840986507265235, 3.4542050240126985, 2.7562255272510354, 2.8458869784303977, 3.2189187383340037, 3.3430886092266623, 2.2094819573395825, 1.961959149916079, 0.8771961221965081, 0.8521879963488765, 0.8640378545759512, 1.2576778466867355, 1.47803902092348, 1.0750741545463645, 0.9553545279770413, 1.2256721356431322, 0.6797654498226974, 0.42528322949398134, 0.7132576828679884, 0.05361668191657486, 0.0]),
array([10.348369134537618, 7.86812977586107, 5.429668147034153, 5.786473332152704, 5.789039142222038, 5.356505361423035, 1.6608538372760173, 0.8252780891778263, 0.7709631327655193, 0.3080911997094882, 0.0]),
array([10.030213808221076, 10.042164783262091, 7.918506622524379, 7.032923320630237, 4.804757598560014]),
array([10.845031781076047, 6.64701577965455]),
array([10.566139074064269, 10.611313578869918]),
array([9.404726985504055, 7.810214645272196, 10.289250771463422, 8.121570730082798, 7.16834113249629]),
array([7.544615630955075, 8.27501824060769, 8.700775774873259, 10.18154943274531, 10.201271003975346, 9.819955472826472, 10.47110128875865, 5.840825733939187, 5.98516623791064, 5.8524098552683395, 7.213078659694377, 6.242401610871511, 3.229792837154727, 3.548785544319777, 2.903832967719721, 2.1853218603240796, 1.2128635854317846, 1.026727067089697, 0.9994758415199186, 1.0622612635175501, 1.3173746784948883, 1.4836562246854526, 1.598903808519679, 0.4147801473807277, 0.2592788934515782, 0.6131110740726822, 0.6783901232714734, 0.40556561539884617, 0.6929519000281905, 0.0753351214457779, 0.0]),
array([9.819665639980341, 8.944095553016423, 9.797221219701314]),
array([8.221597375651623, 7.493483536523006, 10.226914701396874, 9.64728481700607, 7.008767046455144, 6.672534630471262, 5.360177894155364, 7.138372217622591, 5.854274679101882, 4.428473440496303]),
array([9.806270093101963, 10.063288981280968]),
array([9.833643716284435, 9.872954398776583]),
array([9.822873552835071, 9.452309918069878, 9.422441518376223, 9.654585967712263, 8.051258436989555, 7.83708845768963, 7.827693827425771, 8.452143354281397, 7.591686913249186, 7.953697371098394, 6.411486198210603, 6.019034338824996, 5.646801518855954, 5.88850969056322, 5.7314073538211066, 6.052630298517973, 6.62479675195424, 6.170058786301561, 6.386287695295107]),
array([2.8364762451586105, 1.1304848265933372, 1.27334155609289, 1.328137220751001]),
array([5.358363143028766, 7.212829020323636, 1.7962342906273352]),
array([7.907368082313221]),
array([9.305994074235239, 8.778973404757803, 6.562178043791535, 5.872703599739234, 6.830050798974469]),
array([9.104412858982744, 9.489849159238984]),
array([9.393475980690006, 5.660770068246049, 7.160404426964907, 5.436328608707989, 4.893026811817899, 3.96894112887552, 3.9437358566455494, 3.2371505690699554, 2.912310350566231, 2.712878735899275, 2.7800293844169737, 2.9130701126570187, 1.7333960288432322, 1.4377279981381346]),
array([8.981387011778473, 9.282122825593277]),
array([8.999276625878442, 6.957974128679959]),
array([7.054310475978276, 4.062556616327503]),
array([6.979905746139259, 4.6431342906471755, 4.859608657449135, 2.6007782028524966, 2.8413758672281153, 2.623222827836625]),
array([8.263313287166056, 8.167383247519119, 8.438758220238151, 7.494294531905295, 8.405142761958889, 6.9451574937020455, 6.520169167602458, 6.47909239229769, 6.378475175799245, 5.861785269676924, 5.410168409211959, 6.014647740922156, 5.330636138683001, 6.209304074831006, 4.933389446921274, 3.7619723671479024, 3.726203728891515, 4.027298158629839, 2.618388412482135, 3.5105729607361607, 2.91837857288568, 2.944792247320164, 2.9679998527468863, 2.8721590162813597, 3.37694913226107, 2.7986745966466575, 2.798728447839504, 3.129595763618845]),
array([8.364213656681255, 5.517302044849768, 5.760678736804368, 5.70487041981068, 5.865845092670548]),
array([7.546863288687954, 6.637502213663001, 1.4721264411114428, 0.7800856908020797, 0.0]),
array([3.1998788400159954, 0.5171767397157563]),
array([7.446803855669864, 8.110947468284767, 5.655206019688942, 6.036604895170672, 6.237568724531061, 5.807544049882839, 5.000925591611181, 3.637686085780425]),
array([8.327687556897144, 7.116001808216891, 5.645682412208108, 6.212669772675671, 6.13321791900699, 6.700014316803185, 5.774231514877689, 5.721883350819445, 3.6686823431821027, 4.438160889635174, 4.25894157976868, 3.2553369094129168]),
array([7.355232601481664, 8.09551477871356, 6.502181498704923]),
array([8.426870545120448]),
array([5.355543452991255, 7.144552235123289, 3.1583997926641416, 1.9610252783836213, 1.973462384536956]),
array([7.674188645922638, 7.976596802192302, 7.95638194528268]),
array([6.034315148868686]),
array([7.376527968514688, 7.984799366432579, 5.737206726112646, 6.225956773624099, 6.981083197969829, 7.131408498303953, 5.776689716946329, 7.162067171256032, 5.8356154521749835, 5.473378047777095, 6.114492572438456, 4.557957346437458, 3.756942192857961, 4.244616848730402, 3.6385791674300645, 4.0753420047452735]),
array([7.289460668093944]),
array([7.218298844364307, 7.055153996147921, 6.999560769011954]),
array([7.304376975584782, 6.058387481739697, 1.8740617895496552, 1.9164494083765682, 1.2325244079267272, 1.5375938709076347, 1.0898613886018715, 0.12242875907603853, 0.0]),
array([5.347778141281529, 6.996607244764387, 5.657947204144692, 7.064492106759471, 6.279630554749273, 5.74867342495037, 4.9254288059416105]),
array([6.92159969868748, 6.883794262516359, 6.43594043463363, 7.2393350343966825, 6.861460992810473]),
array([6.08276772400166, 5.709389095599098, 5.524737564314863, 3.869343486574402, 5.015694868189932, 2.5063642307916876, 2.40899350528963]),
array([5.553186658095154, 4.588607739088932, 4.774303741950236, 4.379278938931106]),
array([5.336129545392529, 6.165746610494497]),
array([6.080540188845161, 6.870115850500762, 6.9903732210290705, 6.5407519713218285, 6.673169261764429, 5.692876467636928, 6.015034034501116, 6.0224988209296635, 5.5003547256287835, 5.562487047181107, 6.775853896876815, 5.916383939512665, 5.82247779930681, 6.650167363544718, 5.656465904232658, 4.253521162960882, 4.9580557500625, 3.740440546419687, 4.977453344714322, 4.8058735912251, 3.7689773134684432, 4.468700319149095, 4.4503160449007915, 4.993642123888422, 4.255288002564275, 3.253979690228862]),
array([6.133235769137126, 5.618757517734862, 5.288616307413827, 2.7763866942184685]),
array([2.9814715013520736]),
array([1.5342271005674246, 0.6398371625179619, 0.43997176669014476, 0.0]),
array([4.125088966004539]),
array([5.48893740638427, 5.867109894542504, 5.918985684377939, 6.207802396041795, 6.0062286650150885, 5.334637267845666, 5.916223286168696, 5.674281944081827, 5.0765530864209385]),
array([5.434129645596443, 4.617359439837156]),
array([5.1001441679997885]),
array([5.457520931713489]),
array([5.64455559641679, 4.3082913989516465]),
array([4.834729631720807]),
array([5.584349793605709]),
array([5.202007619860692]),
array([5.011494730361162, 5.138625353575804, 4.908216169584274, 4.93056574909338, 5.020517177504168]),
array([4.138204238157576]),
array([2.2103257493953983]),
array([2.902024698911366]),
array([1.3085599441650386, 0.0]),
array([3.690951176735221, 4.863220874161171, 3.9045666799644936]),
array([4.5062644413505435]),
array([1.8762419228637852, 0.5361720154525956, 0.0]),
array([3.9496116370962127, 4.017655152712317, 4.395259387163779, 3.6210537124064848, 4.289245288209252, 2.8681632982771026, 3.174629636561655, 2.8623110465745123, 3.067705583332776, 3.1378772344264667]),
array([4.33006965276703, 3.759726135964009, 3.599387502534763, 2.989767324544059, 3.0054920017622653, 2.7239243585540134, 3.384208268343764, 1.2102066913134926]),
array([3.154959161000046, 3.3314522002513374, 2.1713593116350345, 1.7614569232975732, 1.7018232696595044, 0.8960294145935873, 0.4805091191265024, 0.6247846920798386, 0.0]),
array([2.9722188262395366, 2.727813286949045]),
array([3.4243739922843903, 2.783090167683357, 0.32526254668408, 0.0]),
array([2.583781336335109, 1.3147866970559392, 1.0688700181990562, 1.4319576214641805, 0.0]),
array([2.1099034991036656, 2.4615547604858024, 2.0586219502611316, 1.2021734691202055, 1.50005189625124, 0.8575131443712218, 0.3637262183452509, 0.4924304799832799, 0.6079809834133219, 0.7485551324433478, 0.2783128449208002, 0.5094026963867664, 0.21671454296134185, 0.10505893434063074, 0.0]),
array([3.653377032261239, 3.9297050789673618, 3.550231464318761, 3.5575594509955937]),
array([3.7576198730235535, 3.3821295238699824, 3.318591253088622, 3.351386942432461, 3.453697491045725, 3.135970206329586, 3.3478870647876047, 2.0300674841695985, 2.483528666365374, 2.117199762826595, 1.1508563304106016, 1.277094057608402, 1.507985099229658, 0.9701572356332893, 0.8807127396186292, 0.7592151556385612]),
array([3.308573927978916, 3.3931072932532316]),
array([2.2155548955716013, 1.3944781919870226, 0.8499846022464499, 0.21806037330254013, 0.1835895005773095, 0.705189958054458, 0.0]),
array([2.4076180718882374, 1.7489611854267237, 0.0]),
array([0.8366363826966744, 0.0]),
array([3.560505541321689]),
array([3.6098137641550547, 2.6439247625378806, 3.4205979500765658, 2.2061821256468286, 2.4915730170923927, 2.573858906916524, 1.9508669285052438]),
array([1.4346156855910137]),
array([1.3346254037114655, 1.7618604817348544, 0.48755563440042626, 0.0]),
array([3.502519541285379, 0.9711466377373901, 0.8909441914124165]),
array([3.3165842422992386]),
array([3.445863058788405, 0.7564612884274546, 0.0]),
array([3.378291964067331, 2.6440104078538478, 2.6359738055299484, 3.0404823475154847, 2.575736927674048]),
array([3.3893710597203546, 3.223552077773998, 2.351725353043831, 1.8926527276010612, 1.905055173351455, 1.3929517469496657, 1.776491291459605, 0.7268395415882807, 0.49876430332300253, 0.16841723872921122, 0.0]),
array([2.935843244559597, 2.6835586287645588, 1.7268002628587007, 1.4085728239695199, 0.28818508758172434, 0.0]),
array([1.64147617404777, 0.5368048015586345, 0.0]),
array([0.9496878114980591, 1.4044277401834107, 1.721081819898308, 1.7935814847800815, 0.5621670193654019, 0.6209033942632591]),
array([2.83959557532402, 2.8718007050219, 2.709752809123863, 2.759434180091272, 1.0743489114696945, 1.2653733672782086, 1.2479685713263815, 0.8221162435414671, 1.5733254195654172, 1.4431781591273665, 0.18568384878680777, 0.41585824032760227, 0.6106761181325557, 0.47286269570622924, 0.0]),
array([2.9154159132782524, 2.9126882877666613, 2.8207173308793467, 2.506589709938164, 0.9895325323097207, 0.7883264131947167, 1.2221693104807745, 0.9695298095783129, 1.4387802619787018, 1.4211227594568068, 1.0458233694790424, 0.9182274695878326, 1.7342536929196384, 1.5500510669787475, 1.3843301502944434, 0.7375990479024843, 0.34761344543701533, 0.6663458495764063, 0.6089280806224473, 0.13847111514657195, 0.74811253937062, 0.0]),
array([2.8008378166708394]),
array([0.5448732840896738, 0.0]),
array([0.3901547446359509, 0.0]),
array([1.9466137490044773]),
array([0.28364111453331875, 0.0]),
array([2.1017138710105416, 0.9312300837241415, 0.4586885971173392, 0.0]),
array([1.0875734407468627, 0.3516002960519643, 0.13056996623058859, 0.30909274653929336, 0.20162853913039835, 0.08441433747596158, 0.0]),
array([1.4435607877975134, 0.2271226906774939, 0.0]),
array([1.9982666532672972, 1.4574484770986198, 1.4937838030833803, 0.9708847751843908, 0.24279293663904866, 0.6923749271604879, 0.6706474831118249, 0.0]),
array([2.391202959722556]),
array([2.5889813467524077, 2.2021723701083777, 2.573871933365165]),
array([1.1849538312217551, 0.0]),
array([2.4301181137280765, 1.7536264852724812, 0.8424280704136821, 0.9365564904581032, 0.7491872638930154]),
array([2.1733582291716673]),
array([1.9862482942056094, 1.2023358248712772, 1.1912901846650976, 0.4879252615911172, 0.0]),
array([1.8636104161280163, 1.6218079049841756, 1.413991982003791, 0.9204521780267853, 1.5546741226715297, 0.2675243698606341, 0.3890473994116744, 0.45316806502555296, 0.12494776541242604, 0.0]),
array([1.4144003545040515, 1.028549808613612, 1.6082028624828146, 1.1712663453623984, 0.7632071603762662]),
array([0.8293214573743337, 1.2660428251572586, 1.2260192884021568, 0.5617645025729769, 0.18968013284486795, 0.6774056701989404, 0.23996497455046906, 0.4419424125061152, 0.012035656151036384, 0.0]),
array([0.973034158800958, 1.7993812456743983, 1.2367277640580938]),
array([1.1396337196875947, 1.6258187045078996, 0.2822219920950969, 0.3862641645478956, 0.0]),
array([1.5846488642166707, 1.5421904622008504, 1.5278347251176885, 1.409563615287685, 1.4008211008729898, 0.9302956526860172, 1.5066837155969641, 0.7901169048974624, 0.5037280755906308, 0.44864297777606416, 0.44997146662098536, 0.032486142418189645, 0.0]),
array([1.367652574007642, 1.3565800099508678, 1.2747999908332963, 0.45803984531417535, 0.0]),
array([0.07263009574647536, 0.0]),
array([1.5164288546632807]),
array([0.694601318909615, 0.3231251054807914, 0.28671001032313553, 0.5214276056816999, 0.14426745405508024, 0.18764301444953924, 0.0]),
array([0.9299196852071072, 1.0293375992037548, 0.8974047787046666, 0.9769032747883071, 0.853945646366881, 1.1877379116458742, 0.7269421437410625, 0.5045440612212788, 0.7465582907805799, 0.10514977264400513, 0.0]),
array([0.8682430354143598, 1.0101462063197966, 0.3596347288544015, 0.7569430072169195, 0.25091622784865586, 0.5901838088067355, 0.5770012417260811, 0.0]),
array([1.1527424982941328, 0.4083285591298196, 0.5123180999992064, 0.5130571232062202, 0.737767908291109, 0.0]),
array([1.1900603690301457]),
array([0.5072849505809801, 0.22405166419933398, 0.0]),
array([1.0797301949547902, 0.8487967926141183]),
array([0.9080043753580451, 1.0436920883157488, 0.3085891140617552, 0.5686653046597057, 0.22266554508367653, 0.7243849338014462, 0.6112962135021616, 0.7219984092110497, 0.0]),
array([1.0598170549220296, 0.6112416614301223, 0.3176976614851973, 0.36831474619753096, 0.7257183298526073, 0.5513032167790695, 0.10056334308442708, 0.002450757822248517, 0.11790653612596569, 0.0]),
array([0.9255108744083702, 0.8755311578627859, 0.0]),
array([0.4991284199644279]),
array([0.3677033083863946, 0.5353635146121435, 0.6237041265467361, 0.7391376306652646, 0.0]),
array([0.2672649155934014, 0.1723503391039859, 0.0]),
array([0.4561944831360353, 0.572531954898781, 0.4027028180396897, 0.379073731985245, 0.0]),
array([0.3582450867588561, 0.272410250233396, 0.0]),
array([0.3179794857168655, 0.0]),
array([0.5018557680157781, 0.0]),
array([0.3972758932403577, 0.39348150836957985, 0.4043133436594251, 0.44517692612969323, 0.0]),
array([0.3506794429469842, 0.38619803619790305, 0.0]),
array([0.19610555695096663, 0.0]),
array([0.011643002737622195, 0.0])
]
d = [data_1]
names = ["77"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T8', 'T9', 'T10', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T30', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T40', 'T41', 'T42', 'T43', 'T44', 'T46', 'T51', 'T52', 'T53', 'T54', 'T55', 'T57', 'T59', 'T60', 'T61', 'T62', 'T63', 'T65', 'T66', 'T67', 'T69', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T84', 'T85', 'T86', 'T87', 'T89', 'T90', 'T91', 'T92', 'T93', 'T94', 'T95', 'T96', 'T97', 'T98', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'T109', 'T110', 'T112', 'T113', 'T118', 'T119', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T130', 'T131', 'T132', 'T134', 'T135', 'T136', 'T139', 'T140', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T151', 'T153', 'T154', 'T155', 'T156', 'T157', 'T159', 'T160', 'T162', 'T163', 'T164', 'T165', 'T167', 'T168', 'T170', 'T171', 'T172', 'T173', 'T174', 'T179', 'T181', 'T182', 'T183', 'T185', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T193', 'T194', 'T195', 'T196', 'T197', 'T198', 'T199', 'T201', 'T203', 'T204', 'T206', 'T207', 'T208', 'T211', 'T212', 'T213', 'T214', 'T216', 'T217', 'T219', 'T220', 'T221', 'T222', 'T223', 'T224', 'T225', 'T227', 'T228', 'T231', 'T232', 'T234', 'T237', 'T238', 'T240', 'T242', 'T243', 'T244', 'T245', 'T247', 'T249', 'T250', 'T252', 'T254', 'T255', 'T258', 'T259', 'T260', 'T261', 'T262', 'T263', 'T264', 'T267', 'T271', 'T272', 'T273', 'T274', 'T275', 'T277', 'T287']
def get_taxa_names(): return taxa_names