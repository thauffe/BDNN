#!/usr/bin/env python
from numpy import *
data_1 = [
array([34.72013980205338, 29.890408683995453, 29.499197414842218, 30.622875939226184, 30.220892438930385, 33.196189839851705, 34.25185443984271, 34.41523560330882, 31.102560165068706, 34.177122519831045, 29.646389236331437, 33.57951811683799, 30.661220570019243, 31.455950841094197, 32.99967433589793, 31.33886310422575, 34.89599010915032, 30.566146658774567, 34.002558956440105, 31.448042264320755, 34.31788061872041, 31.69170268165514, 30.744795642842856, 30.37982213672417, 34.20969765562346, 32.23431822571344]),
array([31.74941159393969, 31.79768222435559, 30.054113734668558, 32.65163217350624, 30.333471923319333, 30.43738883572309, 31.703370960616937]),
array([32.36003732051772, 30.124342016519158, 29.452309820290647, 32.207095458493626, 32.43680688290364, 32.07021933401167, 31.15269673839322, 31.085729129574148, 29.91768443256907, 32.25292619138516, 29.355732337476685, 30.84452197662808, 31.799157940340393]),
array([31.373752368608926, 29.487182110935684, 30.77366988097473, 29.970868561388663, 28.231391392544943, 30.771902717273925, 30.06395145606764, 28.623465131337227, 31.35882783472561, 29.201854111719026, 29.89019252900664, 30.3015887123262, 29.775069450092307, 30.00240686990453, 29.58213347198761, 31.08785384164752, 30.260646770531952, 28.33495053805413, 29.761997612564567, 27.88350606621205]),
array([29.359120602822337, 29.692268611195388]),
array([28.168341841973703, 29.98207773692795, 29.483277202516707, 28.78115668087772, 28.557034222750023, 29.382893655471747, 28.6521622674198, 29.43619994298722, 30.022113404977798, 30.19384971067271, 27.034904765809284, 25.93767563535512]),
array([29.80804969738815, 28.12572390236003, 29.70784371756189, 30.07407518336973, 28.95879308436681, 29.130228394449848, 29.49715838105984, 26.39454048512539, 27.273130639556754, 26.504296310166467, 26.58156034449455, 26.98035806330515]),
array([28.183270206127983, 28.45783104429646, 28.8029421192778, 29.271585958076706, 28.47111930115523, 28.49246228346523, 28.894945517293724, 28.13321306831322, 28.194586764845404, 28.69696213223033, 27.28256436931361, 27.485568572377968, 26.246540181573295, 26.20047480262878]),
array([28.139889218658727, 28.138474985983077, 28.26400707837666, 28.073749759076115, 27.11369393315241, 27.83035544812309, 27.71776008883806]),
array([24.922331698123184]),
array([24.728366418536925, 24.14550314805059, 24.12733291521756, 25.670158917405406, 25.32940938523973]),
array([23.921883044185403, 23.74022709031989, 24.90522938845882, 25.076814226328437, 23.554109438647643, 25.819904185820814, 23.389933883846588, 24.70144018607597, 24.427380389061966, 23.20650046118506, 22.74202639331268, 22.678317512235896]),
array([24.058938823869873, 24.210551459575093, 23.701538112064366, 24.34940663407104, 23.948474322464932, 21.103886138323123, 22.351009874927936, 21.077665043748123, 21.858848834137493, 22.17261324853766, 20.769315880325717, 20.94143771989112, 22.295019703633503, 22.31469158769397, 21.8766566301637, 22.73185331701632, 22.059966157134188, 21.75163843009095, 21.27425026846555, 22.855282949831818, 22.614476364331452, 22.04755250164323, 21.814580786990504, 16.007594954176362, 16.075956709972754, 19.900350758759696, 16.521420401689486, 13.926713546972834]),
array([23.568649632102925, 24.80388634977087, 23.66544521632519, 21.400635025675886, 21.5147192539294, 21.95957913789823, 21.656329819826915, 23.000309776942935, 22.33705969931945]),
array([22.438580594620095]),
array([21.383346694666592, 22.597853482394523, 21.282384923417037, 21.900055673884903]),
array([22.21242205921535, 22.270489373061775]),
array([21.517368648672548, 21.334952137115017, 21.550823525876883, 21.31699190136895]),
array([20.508827906893114, 20.871968144300745, 21.191476168177417, 20.690367652957054, 18.72988132483319, 20.00593250872414, 18.98914257242433]),
array([20.556435901657746]),
array([20.48043811718581, 18.174022099216614]),
array([19.153531300458006]),
array([15.14454725660352, 15.395026210843486, 15.467322756016346, 15.250218406155321, 11.655059623076138, 13.760182548496706, 13.266576263750865, 13.35042801967959, 13.551514004918555, 13.123660111256775, 13.09558980605014, 13.683607077222415]),
array([18.894048228470762]),
array([18.073403735698516]),
array([15.418401111132038, 15.749968219595667]),
array([15.5527360218621, 14.100864738563695, 14.598333268247226, 15.4158494151291, 12.38180331932871, 13.814543032758609, 12.354090631610147, 11.72500250765435, 11.842116765688822, 11.722137046022734, 12.91619471463225, 12.719186815698322, 11.026795068139648, 10.78433822652385, 11.456249067826906]),
array([14.041342700113018, 14.454005499122273, 12.155967821966676, 13.093193967508295, 12.34946799502101, 12.488720552242198, 10.549957643435162, 10.555636923412354, 11.31944963050621, 11.046441626984857]),
array([15.985489579883726, 17.513468102175533]),
array([16.08701471953952, 15.597948941040107]),
array([14.369348724014053, 15.543291576037985, 15.264147309295236]),
array([14.968617796096268, 14.181581273422092]),
array([14.79466652383997, 14.361049264549456, 15.715102338234756, 14.104953182845579, 13.742602510053201, 13.757970435374672]),
array([13.881319888573119, 14.060331926062236, 12.125360585664602, 13.10816646822199, 12.266630136836858, 13.601110104194825, 13.298151731133142, 13.706600922747528, 12.391452907085952, 13.485883210470348, 12.861637479322047]),
array([15.412472749853444, 13.944557257280008, 13.280071213537399, 13.23233565838303, 13.724523898554677, 13.503253514188525, 12.851002968536358, 13.3404964437187, 12.75432311270387, 13.63583706851784, 12.829131290453441, 12.827461539177627, 13.158848540340784, 13.480599237902634]),
array([15.092978810859885]),
array([13.20791203443425, 13.090135964886192, 12.975079053263558]),
array([15.504596976815249]),
array([15.494776524441022, 12.494930239765653]),
array([11.3134748301372, 11.584025780353443]),
array([14.066348065975218]),
array([14.637869830311608, 13.941296013920983, 13.424072675997179, 12.789836948124035, 13.701985979645883, 12.650401882110945]),
array([15.130386716250465, 14.14930052926869]),
array([13.135840086966521]),
array([13.433193345555335, 13.55311678026249, 13.797131969549794, 13.66562993683733]),
array([13.666441252720336]),
array([7.767091042406024, 10.721641545470082, 9.97484721810233, 7.6509995745403, 7.516752393166989]),
array([11.695992061662752, 12.011496807374048, 11.636079816418746, 12.839853591334085, 12.443246538725754, 9.948817834149914, 10.179688650291949, 10.478342970340902]),
array([13.010640746128198, 13.118367096166002, 13.59000020299046, 13.201649187764861, 12.692372948992913, 13.218891721105148, 10.873193578236343]),
array([11.75047388854123, 13.122745959249617, 11.817857791607876, 11.110056741590409]),
array([12.911622656195002]),
array([11.907336491136919, 12.330264436150316, 10.828321668553935, 11.270656116309542, 10.91461956940093]),
array([12.332757306705695, 12.525762451513287, 12.587019650155378, 12.0643240284583, 10.138437142959344, 11.409444278161077, 9.538539460779454, 8.093970738004442, 11.14489056942105, 11.011960717163266, 4.972080988832833]),
array([11.984245962999694, 12.291146788026131, 9.24945085884424, 9.590999716561683, 10.732921553228055, 9.952169142587096, 11.29935047406542]),
array([12.541868526304711, 11.422878215024934]),
array([12.731827289888619, 10.916077922885835]),
array([12.13902545364502, 12.623316783225812, 12.678954411770492, 12.072364619804626, 12.637894838027739, 12.235068670789381, 12.05916105574793, 12.45480675120095]),
array([12.047857791365942, 12.400008102308547]),
array([12.536874725176805, 12.324672553552476, 12.344851042109147, 7.565527502128829, 10.984206067689398, 10.568164644453812, 10.37749895418888, 8.365825192304568, 8.342451545145153, 8.484969157215962, 11.446501388544092, 10.238152365534042, 10.459538680396118, 9.891621373523929, 11.477552879077887, 9.513477770654488, 10.62651372161566, 7.508515903382332, 11.430238646937326, 7.321956467228199]),
array([12.28987276828917, 11.987141294786772, 12.514493713459382, 11.802776175233983, 11.653094905292187, 12.549188291718059]),
array([9.6762084739813, 8.83020908283747, 7.6343178538113765, 8.788096033024221, 10.325184920240211, 7.351862938754212, 8.974124592965786, 8.492826945386065, 8.039402276014219, 9.016392303887656, 6.662324380461882, 6.491816508157, 6.54290278472721, 5.4717876123981135, 5.4901339006471055, 6.113130877961843, 5.616587163697847, 5.383359123446288, 5.295193566489639]),
array([11.770883889617878, 11.827040026695194, 11.747063352209324, 9.475589921589078, 7.72714528285289, 10.824397623288991, 9.677274105290893, 9.764757246886457, 9.563943121767103, 10.273083296724916, 10.562079701919412, 10.391414446441457, 9.290805253452408, 8.11785658663118, 8.669555968150425, 5.52961003235505, 7.0678252825086485, 5.876485393700216, 6.739937924972115, 4.236629787186008, 5.053895255015134, 3.7248030569361026, 3.123198528181835, 2.7882726452555175, 0.8164817396175473, 1.0231739594116647, 1.456248530309315, 1.087754551890078, 0.9746340302670238, 0.24532391926553054, 0.6924592731618789, 0.4779274826798575, 0.6440911476236133, 0.750932968488605, 0.0]),
array([12.144042057514568]),
array([9.151024684543454, 7.534780194534336, 9.65187126511657, 7.790028663780234, 9.967572428171334, 8.685533412197987, 11.289716620010077, 8.847182244181248, 8.77838883618404, 7.146760006623245, 7.104387622070899, 6.455081985009953]),
array([11.685184964568592, 8.551203281695756, 7.345389838051646, 7.8556517107024835, 11.12780612786302]),
array([7.345164099268934, 6.553083515779573, 6.340131150770852]),
array([9.194974651799152, 7.949974893235516, 9.233988526730036, 7.481501835798179, 7.806914705999448, 7.274399354432445, 3.9619984062472398, 3.839043782980255, 2.7266919434392545, 2.9940594324599004, 3.323756939590662, 2.42841419132109, 1.2003754693319046, 1.6202270418599691, 1.519628213228006, 0.9693184949666069, 1.701090166080197, 0.4441176245407077, 0.743800792728553, 0.0]),
array([8.762086213884404, 10.225814098100143, 7.310026509175029, 9.742801667308074, 11.03341511295488, 9.03382212486893, 6.712197009521829, 6.30241797480187, 6.010283343927058, 5.5364384095129395, 6.217437785411609, 4.939626209567193, 3.1277604121005385, 2.450667343007983, 0.9846767111046695, 1.6293652243386312, 1.5163029108709067, 0.8147335897498014, 0.6950347782242214, 0.0]),
array([9.611018337993594, 10.354085519020414, 10.372491988112456, 9.904985139290213]),
array([10.271210785518361]),
array([10.106542686564373, 7.309898644060002, 9.18734742116282]),
array([9.839153643436289]),
array([9.558982540236537]),
array([9.101272693013982, 6.795555554960693]),
array([9.032881085341593, 8.426364071508239, 6.515866462949445, 6.249934121639144, 6.500816781133536, 6.120827300313675, 7.16002898116197, 6.526962166169468, 4.43930936108787, 5.0590025558014515, 2.209918248870873, 2.043677861452948, 1.877138824494186, 0.9789058319782548, 1.5465987150038092, 0.9637551058295509, 0.4548050219719514, 0.6532447602818843, 0.0]),
array([9.027954503174497, 8.176344934917566, 6.790770262737402, 6.8683953099334625]),
array([9.521856768107131]),
array([7.5582777490435005, 7.808479817783391, 9.017143221144282, 9.050132858834756, 7.8747810399321265, 7.914293387678435]),
array([9.252287326448055, 8.31081775947208, 7.657564359176633, 7.561552134505863, 7.788415692949637, 9.079978947549707, 8.88832401446046, 8.796075816425605, 5.757710330232686, 6.535718642508763, 6.156338797780868, 6.248952040171015, 7.02450138745479, 4.860618074529003]),
array([8.814078276484945, 8.973461767777172, 7.394486380860705, 9.278639786815182, 7.605860917765722, 5.542844427177862, 5.012966508337105, 4.081988739968418, 5.025493531021771, 4.572028002633232, 3.2523018470844867, 3.5984792715814797, 3.3777693438141694, 2.9721816592907033, 3.224798826970634, 2.0016741905828463, 1.9888028026668443, 1.1831572772042902, 0.479226302534629, 0.4217873656540533, 0.0]),
array([9.149404509269232, 7.498427910295348, 7.7397972071641, 8.122375736365491, 7.859252682712089, 8.383449192374206, 8.717254295310681, 4.2159737997401745, 4.800288601857064, 3.199667301590241, 2.8439404316666943, 2.4309232257818976, 1.2940213961941973, 1.1248131785295525, 1.0042315063673097, 1.2313640415197726, 1.1955478310921412, 0.9472152734876845, 1.5361039962826857, 0.5399011131640719, 0.752514950403343, 0.0]),
array([8.99148379527815]),
array([8.598992698991797, 8.002747092193564]),
array([7.383344492402663, 7.369971166778252, 8.565480131734516, 8.178139947520773, 8.46756857861888, 7.34966499572603, 6.1729585931926705, 6.117947299000197, 7.037963142170291, 5.816025152541826, 5.956143996010857, 6.583918677723496, 4.661219713471528, 4.929389081393063, 3.9051749264953175, 4.520704838208597, 3.5405667469575244]),
array([8.414196682639751]),
array([7.304448117098208, 8.229703470884862, 8.043047641681921, 6.114779236904147, 4.018863697775959, 3.82272697427611, 2.9507142555451686, 3.073370101415354, 3.3028932765712877, 2.6846704350585595, 3.0500139282072714]),
array([7.699491341081143]),
array([8.224759849465977]),
array([7.24731497409428, 8.090746929434879, 7.2054297026746355]),
array([7.710350014191061, 7.641862618905344, 7.802986695334495, 6.6885341175750845, 7.136758111618131]),
array([7.15904168858521, 6.500754974563895, 5.723258567927948, 4.6900259105991475, 3.25654985670318]),
array([6.853792599798622, 7.187170283843956, 6.000548738430927, 6.527461734465103, 7.146742848015047, 5.940031553948522, 6.424588173469642, 4.520699005156736, 3.6054320667136333, 3.3888680927562307, 2.3171137009269223, 2.3464759585669697, 1.3559995156826403, 1.513384860907836, 0.899179506335297, 1.363409101918015, 0.0]),
array([7.73488800093631]),
array([4.361710665008986]),
array([6.863045697365259, 1.3356810911113626, 0.0]),
array([2.466555174398023]),
array([7.575599557960353]),
array([6.93418791066248, 6.2872810855256045]),
array([7.227275097355027, 7.206262164762185, 6.881965513189504, 5.978218067906352, 5.35142274264379, 5.074811187130304, 4.454153539644726, 5.103866891292295, 4.137972199659609, 2.818640415131128, 3.223509473972083, 1.5513573112420782, 1.3017954722858316, 0.8183769920936099, 1.4927010853558484, 0.7786122381488928, 0.3506925422902483, 0.0]),
array([6.84447387255185, 6.307453647472092, 6.507000918883368, 5.897677095650359, 6.661169465430125]),
array([6.706691808657622, 6.590382607207779, 5.471507056327045]),
array([5.785849448117762, 6.343010869375823, 6.241183131183681, 5.839564089611015, 6.399289654599788, 5.848168728379722, 6.899514444017241, 6.40283928574393, 5.5500582755967445, 6.7585166399130845, 5.753377889550924, 5.98406390459911, 6.883392153081358, 4.948947066498, 5.090448184866317]),
array([5.6767774152275745, 4.720913939591636, 5.150110460686529]),
array([5.454074719089911, 6.577739867362743, 6.281124721201294, 5.650065020864307, 3.839399819447482, 4.7744890514452765, 4.42349452800534, 4.87351851257654, 2.2300345918044515, 0.975583726762954, 1.2014802737789643, 1.6113559789022476, 0.8161046395879565, 0.22769041602513418, 0.5325005239586093, 0.0]),
array([5.365483524806204, 5.853979932873635, 5.434697659812258]),
array([5.793641781876464, 5.186508533140264, 4.808783738661944, 4.422453080840585, 3.065752329581298, 2.1552150150407936, 2.476973874133567, 2.501784170326085, 1.9065302941989803, 1.2267992686063534, 0.7972790335952791, 0.8751564667207883, 1.2405573742576261, 0.47525061808608765, 0.7039942766198467, 0.004399985700793474, 0.0]),
array([6.162460543422166, 4.063633918995505, 3.4068000275353603, 2.759732024227066, 2.8830068173021255, 2.554502569568099]),
array([6.090789809573866, 4.834997033565499, 4.55380038950015, 5.075734517042451, 5.19602581554234, 5.000324655240988, 3.4296661127294072, 2.843321158819575, 2.875471011244223, 2.4406248207824848]),
array([5.37374746974345, 6.0143677574469585, 3.890476688124088, 4.785209329529469, 3.410655496388789, 3.2047768339864646, 3.0103792185104163, 2.7060616361926835, 1.4704408514593228, 1.0571087197590499, 1.7987155738218787, 1.5852386430817058, 1.2471117909583755]),
array([5.456690998507164, 5.370903009403676, 3.7165197785678266, 4.427733277067435, 3.698787933151466, 5.305522574787592, 2.943235456182697, 3.3814931500742875, 2.749160104413596, 3.290099312712169, 2.8584991170542575, 3.1101746207807865, 1.8171509708280758, 2.5215230383349585, 1.5275360646236233, 1.2426755186035785, 1.2908680008821785, 0.5845512074471715, 0.4151888352628541, 0.0]),
array([4.628835312217006, 5.237033413765688, 3.2844392030037546, 2.628506732031105]),
array([3.969628124676299, 4.099053489569851, 4.445310980121634, 4.444668850368263, 3.703119850554298, 5.116179119466049, 4.030274196638937, 2.957616267638303, 3.4751241732641223, 2.2439415911949827, 2.179668376889559, 2.294550348276676, 2.487291298067278, 2.2459850626087157, 0.8855720124235017, 1.2422095042653694, 0.30937491967229513, 0.0]),
array([5.356548849916309, 4.6381167428751064, 2.9955050014686697, 3.217671384252286, 3.536183378974318]),
array([5.437083501708133]),
array([4.124352722680905, 5.069878081604512]),
array([4.983613438924213, 2.440111449144952]),
array([3.613938362499689, 3.7947461150698496, 4.487365925604735, 4.738079724534889]),
array([4.853086142243556, 3.193981027522682, 2.896766001895716, 2.9262797816922586, 3.2448797925380606]),
array([2.170459474640376, 1.8572550903853409, 1.5796841331942553, 1.35981053162619, 1.7836700241105392]),
array([4.6210060903391525, 1.9492362442426392, 1.318099014699695, 1.2035030647951601, 1.1324430068015143, 0.05029289405981541, 0.0]),
array([4.325331400678455, 5.021243481369917, 4.936410261671847, 3.007941080068027]),
array([4.017962276576738, 4.670582974234708, 4.36391389604681, 4.5475518248393145, 3.581734870712677]),
array([4.434577271809811, 3.437444813650563, 2.76399662265126, 1.5664411827205125, 1.4708267255691725, 0.5464692491378352]),
array([4.463138074093797]),
array([4.239401412822878, 4.463598090967258, 3.200629498903574, 3.237095694166467]),
array([3.455034625913186]),
array([4.115414178133763, 2.8489205664338106, 3.4047294813608944]),
array([3.3949300999500482, 3.559635802477044, 2.109063276887579, 1.6446281601868529, 1.3373792989463062, 0.0]),
array([3.1074637517036723, 3.427982173247992, 3.5339591934755004, 2.661752781931919, 2.9458671962237113, 3.4092500769431773, 1.5731710871819988, 1.711262973739903, 1.0627553147486313, 1.190174474459465, 1.331250932501355, 1.4602684688246304, 1.1789528242397045, 0.28399329275952606, 0.1826774346958403, 0.0]),
array([2.7111431570725886]),
array([3.3120398198079735, 1.6614068371992774]),
array([3.0820452833678313, 2.5831201961707118, 2.7705784248771836, 3.0941303040758887, 3.0256305447891156, 2.827242647831193, 2.778128932453568, 3.4646896066482977, 1.9667442490945137, 1.835557130198978, 1.7801579770408633, 0.8517947209684397, 0.9787191509348825, 1.4704058036056469, 1.4908042425249384, 1.2065198631417393, 1.0002296976631073]),
array([3.6570901507561553, 3.2980345004260054, 2.804941539576118, 3.4638824791144436]),
array([3.411538808064729]),
array([2.77111405453767, 3.274552897508639, 2.8750749375827263, 3.299860104377851, 1.9393503118802964, 2.44093844394526, 1.9361718185892922, 2.265983615613239, 2.245674680109171, 0.9438740006979527, 1.1120482582391, 1.0275657589987568, 1.2735464730465882, 0.8129472082524489, 1.3741603958460575, 1.3581736485561926, 0.14091590744939197]),
array([2.8892795927039097, 2.7677317200570624, 2.0954335795857477, 1.244884330768188, 1.0939351374684287, 0.9856469830927932, 0.02201313858105189, 0.0]),
array([3.2679042047482723, 2.843654224405993, 3.094131221334877]),
array([3.2358117404278914, 3.19098353622143, 2.6205757258988176, 2.087598039559276, 2.3192693510223714, 1.8953542209403698, 0.9721198307688967, 0.8038085732992752, 0.7839770148268517, 0.9396541682542292, 1.5353861074568116, 1.7242582488660072]),
array([2.639491834454992, 1.724609660305184, 1.1204527659737649, 0.0]),
array([2.7857101590504927]),
array([0.8833493230854411, 0.8491339220796587, 1.2738202095231763, 0.9909383149727508, 1.4195971998233046, 0.0]),
array([1.1626443217291436, 0.7598880917991586, 0.0]),
array([2.6117526302430667, 2.5873776568499944, 2.452443586311473, 1.9563360480155474, 0.9039144838655354, 1.0038078657956788, 1.4154228245358218, 1.1449109815590597, 0.8974111420772461, 1.0499622476143657, 1.6716598399683378, 0.43556716006287954, 0.5725248404111518, 0.0]),
array([1.2782675992418537, 1.7119946959983563, 0.0]),
array([2.6592429077909308, 2.404902490839662, 2.318229978610976, 2.16897987342297, 2.0044324885459748, 0.9790828034780915, 1.1727447904209432, 0.0]),
array([2.1181552676028623, 2.4438351073766134, 2.0835016993176794]),
array([2.0598962125502673, 1.022395934419177, 0.15888928208312758, 0.4045839547800495, 0.0]),
array([1.6317871884005206, 1.5441766762583635, 0.8047296295409095, 1.2156219119399574, 0.17831500955986168, 0.0]),
array([1.3445367282050567, 1.7317657468734304, 1.3463000226353845, 1.2902047649515898, 0.9355530646520698, 1.5002434447848565, 1.5686520933263675, 1.2524413666315302, 1.0129409131155396, 1.3886259148167177]),
array([1.863607950493823, 1.5860136689900468, 0.0]),
array([1.7495063150145276]),
array([0.9822254536316665, 1.416917626965633, 1.0038546504427681, 1.0498749711604405, 0.851421873089076, 1.6684491938528823, 0.7257480364190156, 0.7744423628833611, 0.12895452405744545, 0.21070187780153626, 0.3069546098429987, 0.15121568701708665, 0.02045758674142978, 0.0]),
array([1.3775159303010955, 1.6884399278752127, 0.836143425803403, 1.690935376739027, 1.2285064751259793]),
array([1.0523271527887887, 1.7433408788273244, 1.6296072108681479, 0.8291940091283531, 0.92567445923236, 0.9648326988035698, 1.0580721142898548, 1.119559035671065, 1.039130360769844, 1.7900886438195525, 0.6327996884567453, 0.7244051472383765, 0.46972052017713595, 0.6384291389657369, 0.3936703196912436, 0.11356322498202104, 0.0]),
array([0.9337789479973814, 0.0]),
array([1.2635183145665032, 1.335268195249319, 1.327347772267525, 1.1016350795294039]),
array([0.9818373401499413]),
array([1.3333405517047523, 1.0749258743955714, 0.0]),
array([0.9642032769612991, 0.9354183847002335, 1.2254151235572348, 0.0]),
array([1.2890489539412315, 1.0751239401036359, 1.2774687284672088, 0.25660059232453725, 0.3508443828024484, 0.5386520051759436, 0.0]),
array([0.9063891776220746, 1.1046180140197546, 0.9435571104720435]),
array([0.9054605819824435, 0.2890996838020717, 0.5426829700032736, 0.29318185342752934, 0.05110669272219137, 0.043511686577230604, 0.0]),
array([0.9724506625238979, 0.0]),
array([0.12108821199397703]),
array([0.6962943838488321, 0.42258905057373286, 0.0]),
array([0.6515937389796748, 0.0]),
array([0.16677941094748805, 0.0]),
array([0.2822849519517162, 0.02549012923876974, 0.0]),
array([0.19582551742150722, 0.0]),
array([0.3444638008559541, 0.0]),
array([0.06883486541614178, 0.0])
]
d = [data_1]
names = ["6"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T8', 'T9', 'T10', 'T13', 'T14', 'T16', 'T17', 'T20', 'T21', 'T22', 'T24', 'T25', 'T26', 'T27', 'T28', 'T31', 'T34', 'T35', 'T37', 'T38', 'T39', 'T40', 'T43', 'T46', 'T48', 'T50', 'T51', 'T53', 'T54', 'T55', 'T56', 'T57', 'T58', 'T59', 'T60', 'T61', 'T63', 'T66', 'T67', 'T69', 'T70', 'T73', 'T75', 'T76', 'T77', 'T78', 'T79', 'T81', 'T82', 'T83', 'T84', 'T85', 'T86', 'T89', 'T90', 'T91', 'T92', 'T93', 'T95', 'T96', 'T97', 'T98', 'T99', 'T100', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T110', 'T111', 'T112', 'T114', 'T116', 'T117', 'T118', 'T119', 'T120', 'T121', 'T122', 'T124', 'T125', 'T127', 'T128', 'T129', 'T130', 'T131', 'T132', 'T134', 'T135', 'T137', 'T138', 'T139', 'T140', 'T141', 'T143', 'T144', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T153', 'T155', 'T156', 'T158', 'T159', 'T160', 'T161', 'T162', 'T165', 'T166', 'T168', 'T169', 'T170', 'T171', 'T174', 'T175', 'T178', 'T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T185', 'T186', 'T188', 'T189', 'T194', 'T195', 'T196', 'T197', 'T198', 'T199', 'T200', 'T201', 'T202', 'T203', 'T204', 'T205', 'T207', 'T208', 'T210', 'T211', 'T212', 'T213', 'T214', 'T215', 'T218', 'T220', 'T221', 'T222', 'T223', 'T225', 'T227', 'T228', 'T229', 'T231', 'T235']
def get_taxa_names(): return taxa_names