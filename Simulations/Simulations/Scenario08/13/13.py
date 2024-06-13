#!/usr/bin/env python
from numpy import *
data_1 = [
array([33.8765858626745, 31.257819593528147]),
array([28.434755588105833, 34.01938875643774, 33.452647951832056, 29.903352235018694, 33.741535237683664, 28.968606855945545, 28.147134661349615, 32.29452611915438, 33.32084476546288, 31.389130421768392, 28.48808079223082, 30.05008986043902, 33.54581905420403, 28.541405118761173, 33.59802538228831, 33.82177506292122, 30.785091174356637, 30.522210815571388, 34.23586107377719, 33.99700505860402, 33.32924132081955, 31.90228633137312, 31.50797802733286, 29.765506013755196, 28.84075581487725, 32.8707957955736, 31.5215888375975, 31.624640905645084, 34.569282770557464, 32.872273277689516, 34.317778292482, 31.190135148511477, 32.25280186534949, 29.302886070519467, 29.584286508045516, 31.09235636314276, 29.221741559405345, 26.332274159230444, 27.57385622598003, 23.25908285810894, 26.027529350636172, 22.926126623466992]),
array([33.47250126381569, 34.44437019760669, 33.45634451804599, 32.55369335700338, 33.51866643008634, 31.825347611028054, 33.32144421240884, 32.21676711526732, 33.93453058069609]),
array([33.4346658343926]),
array([32.69311234906087, 29.790464646941974, 29.46796600068222, 29.23267657851122, 28.104759408880554, 30.830952232363526, 29.82665695051582, 31.15041148900025, 32.16100763679167, 29.931130975747593, 31.992210864754924, 29.805949717768712, 28.520673532256836, 29.83197759761066, 32.760342418016116, 27.270109594509805]),
array([32.74546395570747, 31.278202877196943, 31.697143390014613]),
array([30.56780667092286, 30.917316574849753, 31.18147916070416, 31.51127239523329, 31.14506527603104, 31.618467071579012, 30.601201532916452, 31.62596534978929]),
array([30.230966544559166, 30.379222850730347, 30.17935642968506, 30.63427357097512, 30.46081416972993, 30.89052048478421, 28.41192296519835, 30.764319727386134, 30.69780461416444, 29.282199392166582, 29.982315692655323, 28.843339579205796, 29.535288069396348, 30.900487129047516, 30.602714177615066, 30.518270384359546, 30.582672569150557, 29.674197890030154, 30.307183762435898, 29.339242267876205, 28.9240544922216, 30.43483304573738, 29.47117038459954, 28.414717948325073, 31.039527996291252, 31.135153350035626, 29.920459724545022, 30.367309657057376, 30.47076217389618, 27.66715267948138]),
array([30.552639493139477, 29.688447759480137, 30.988335378580693]),
array([29.710098384698558, 28.80766665064244, 29.8811505972916, 30.548034860273656, 28.81960075917557, 29.685140227233042, 30.275817753180505, 30.388754815595945, 30.256735887970976, 30.076874978804764, 30.646065131400857, 29.515303199939396, 29.204083670603794, 29.145620642207003, 29.263110292710337, 30.688048445137117]),
array([28.807472474284687, 29.60118479774117, 28.27430678775546, 28.90154789003812, 28.395897991988686]),
array([28.64269066067939, 28.34284772291063, 28.163412598965476, 28.616134471178597, 28.63552277413248, 28.417249023243027, 26.550264086523136, 27.068755307694254, 24.936305426517684, 25.19069940694862, 25.538777346400092]),
array([28.446218688720055, 27.581546647459767, 27.81179076323499]),
array([27.405769407128457, 27.37894196886003, 27.732033475384046]),
array([27.369785138189762, 27.19858647565733]),
array([27.260859936260328]),
array([23.23110260964015, 19.08419823600964, 20.237377201606222]),
array([26.600718802036294, 23.22367704186534, 21.87917432318941, 22.11471852192451, 20.12742095842883, 18.343493572604864, 19.101617469671197, 20.25111491586575, 20.21273967081027, 19.58967896326129, 19.529129527710275, 19.24456474556149]),
array([26.786281059536442, 24.81396715967379, 26.531877453017366, 26.03911947330461]),
array([26.505454070770096, 25.230280186725892, 24.01864096609294]),
array([25.43459289493098, 21.853100521470836, 20.958670617507927, 18.120101407192376, 18.841557668520558, 15.650652294578702, 15.05762330116564]),
array([23.843365963319254, 22.98417466288378, 21.15183913043564, 21.819118682587835, 22.602952896400488, 21.92512774658871, 21.669188378614315, 21.87162128304502, 21.723477558929822, 21.840438919859064]),
array([25.518076079797414]),
array([23.827455811953307, 24.38119361498621, 24.453468431971448, 22.189304394846097, 21.644825367216136, 21.395795681860513, 21.369945239459394, 22.562088612279837, 20.265538737827068]),
array([23.244827862403433, 25.38235904792504, 23.878781589876763, 24.031799329718208, 25.610274958908644, 25.456143497620154, 23.888960284113516, 25.169867315674892, 24.591777030042433, 23.242552266241862, 24.632327936677477, 22.91229032853211, 22.982620167438753, 22.795176343308412, 22.637618971286784, 22.45323067963301, 22.57537082075557]),
array([23.04655111311957, 25.170069118487728, 24.030883259129272, 21.94857058297582, 20.864439098591955]),
array([24.932176831373564]),
array([24.133615524426872, 23.02328732121893, 21.1496199211945, 21.06400115678656, 22.44312158196844, 21.612332755221566, 20.505081534686454, 23.015703221521598, 22.91079631786416, 21.209380555093723, 20.51856135819272, 19.05237434898644, 17.960506342890508, 20.055662047868466, 19.310467409101783, 16.637676577229925, 16.039127418870965, 17.904281733428125, 16.84979081647128, 19.08681582871118, 17.560046756744498, 19.99593581479519, 18.774249707566877, 17.616207340876983, 16.48093044454534, 16.562445339672628, 16.471733555960043, 18.991428930219072, 18.13912636717245, 18.382090570617535, 17.24188793904166, 16.24066269791807, 19.39431163676539, 17.77485236337924, 15.888776856660254, 15.13277182688082, 15.19417781335886, 14.573596890567439, 14.128880716065458, 15.020313873624438, 13.991040620008745, 14.325164691170173, 15.486907771477501, 15.263033317902542, 14.38492976455368, 14.9626812910872]),
array([23.614430182590425, 24.111731380389223, 23.53875126027768, 24.35322176616298, 21.05709571990248, 22.218611942750726, 21.592687790929535, 21.931176226445302, 21.187639186071344, 22.010585264857056, 21.563670761189975, 21.184829628283573]),
array([23.30415837915136, 20.834878661723966, 21.774246477127758, 22.79138247487809, 21.47003309297865, 20.603229420966294, 18.866206751130736, 18.05968990723712, 20.366470953624425, 18.75624283694941, 18.02750346268064, 19.6559042567806, 18.429492739205983, 18.04235957987102, 19.692957508956358]),
array([21.31976218926496, 22.96098757723375, 19.978716911267973]),
array([20.034915401476525, 18.702964217906903, 19.26554129771371, 19.379442761125826, 19.656549223484596, 15.228755935889092, 15.489916034210463, 14.891606329663972, 13.069466202705986]),
array([21.10276540205981, 20.72275817051336, 22.54814960920849, 21.46993826311449, 20.231939070791444]),
array([22.90857819812048]),
array([20.83162399549651, 22.926170349660314, 21.640291180900157, 22.2740400540222, 20.440752680302342, 22.978368172835193, 22.369437200050704, 21.743148432105876, 18.84727805929451, 20.22579078295437, 18.400563555286254, 18.54714573009882, 19.280371944854803, 20.249047956377925, 18.247530279209865, 18.467590860916616, 19.718908556187078]),
array([18.787453526713065, 19.37000372180829, 19.81713241172458, 17.79429442195899]),
array([21.73886535666079, 21.284924434214155, 21.714792872105694, 20.820874610868998, 20.547337812259382, 21.786332122450997, 20.59695019389079, 20.793846924142702, 21.444022643174964, 20.82144202469332]),
array([20.51341912511609, 20.571472169432266, 19.362403168949097, 14.599436694853605]),
array([21.121489299537014, 20.979581866361833, 21.206914390403238, 21.436607343250934, 21.142519417632847, 21.137413817591305, 21.309420665026998, 16.528094003434468, 16.884765084572095, 17.459678289036013, 19.756606554291583, 18.515077135161143, 18.07698962570917, 19.929834455172443, 19.344277819499425, 18.243025745641408, 16.353435534930888, 16.17897330409665, 17.278267754365384, 17.4539668436684, 16.021001355689176, 18.21905578143494, 19.86577564269757, 17.68100202276467, 16.040525907413326, 18.441107383848234, 19.6138128126378, 18.388601685930393, 18.674892790205725, 19.219804070284166, 16.395366112989795, 16.236532274127978, 19.425045022433544, 19.316990307077567, 15.500421667324684, 14.620776014067262, 14.680128186653333, 14.022372875823859, 15.34587743220168, 15.758715308612215, 14.952535794610583, 14.147890518550959, 15.365802514082498, 15.296629967171437, 15.548218658801748, 14.103083844131563, 13.921486894759829, 13.441820788612452, 13.814418392200784]),
array([20.957218982138084, 20.282896069533912]),
array([20.63498929654396, 21.038859200832103, 18.821435270599892, 20.374120288945583, 17.26019203687564, 18.92043917704994, 20.07059048773485, 17.834694821606018, 19.418515650972793, 16.65864202159164, 16.548561614722757, 16.301369719355243, 18.41358190149804, 17.203669180801224, 19.54604926269537, 20.205774859202382, 17.12693496577083, 20.14351703691816, 17.469821231647177, 19.982679284235463, 18.136734856833694, 16.149065663564414, 17.164656946569355, 19.37034581330135, 20.06972470332695, 17.09032126276084, 16.098690906142085, 15.856714844560512, 14.921364052827622, 15.110936248782192, 15.641028567544069, 14.723951601995074, 14.224686873040705, 15.924018456961697, 15.478136496134633, 14.583654500458463, 12.774018269253519, 13.340168651279155, 13.138115336554455, 10.323032623746787, 10.492900750884294, 9.510059503618569, 10.514007588662455, 9.81946351114216, 9.52972845433749]),
array([18.929740057021466, 20.224342231532237, 20.042511676773117, 19.24088558957765]),
array([20.812796219285612, 18.124573901676886, 18.471080077400725, 19.554142713493004, 19.242384567547372, 19.02448723351976, 19.428088368081873, 18.22223052626738, 18.149461250854664, 19.449719681651466, 18.249914837611872, 19.355570833794694]),
array([17.859842835316464, 20.306706603712932, 17.332221738745723, 19.643374917019212, 18.55763561029855, 19.888024526794783, 17.60142372040859, 18.441532414950963, 18.53265622707491, 17.23501328736551, 19.6705663532102, 16.41932688029591, 17.3520454810057, 20.18455919704069, 19.281414207906685, 20.038765431199597, 17.679734555934942, 18.57798480478516, 17.35592721936477, 20.05898171126985, 19.43859776610633, 19.49707410571317, 20.189817300393, 19.92567068545329, 17.468180740419292, 19.38401402576814, 19.80023344188403, 17.353612390992797, 19.382159256835767, 16.578387082051627, 16.094115016226535, 20.001302602771787, 16.758504296767853, 16.400529098734893, 18.536225528172235, 17.582387287845172, 16.89953448902969, 16.068170426422313, 16.210842789039074, 16.14040832750551, 19.966191395950325, 20.35926848181045, 16.06936530646839, 16.138215149226546, 17.71925184746297, 17.412589550139376, 19.957382797069084, 18.415618810536003, 18.663466285038268, 19.091921677100107, 17.874518002450383, 16.2418838987135, 15.835825120375073, 15.073637066114973, 15.746977613615687, 15.145770657289692, 14.957529746627422, 15.80405265578119, 15.425946867589788, 15.260520517111953, 15.884774611279113, 15.846755677075235, 15.1327838588053, 15.81148440963986, 15.275854808489344, 15.686365642385685, 15.846147562975787, 15.405239968942887]),
array([19.61528338584368]),
array([19.329075908838753, 19.284696308049842, 19.239094359581763, 18.33155857156639, 17.798072418112405, 18.034118202154772, 18.965970751943594, 19.340601812515224]),
array([18.21467012819379, 15.7534666288125]),
array([17.73262171616257, 18.275295586633305, 17.465884898483846, 18.087136659804663, 17.18798557929051, 17.26605659687229, 17.726867679072807, 17.483263365374658, 18.485268428790263, 16.74161909424481, 18.08869234620187, 19.26715506219434, 19.430425761815847, 19.29076284968651]),
array([18.480121504283478, 17.90190166000801, 18.68516592271076, 18.59060911020728, 16.954976210406574, 17.267186772965065, 16.775676488393174, 17.751626285440814, 14.014562795155882, 15.862091619184422, 15.769799235614018, 14.235423842693717, 13.663083830809871, 13.401984362925473]),
array([18.370334487235688, 15.987835963649196, 18.464963835229256, 18.154730577917192, 16.187955098222144, 14.99665772839049, 15.82091402512977, 14.486495418382706, 12.536929955346276, 12.352240813926734, 13.55601490338626, 12.011237002669862, 9.5270612987272, 7.300374596354842, 10.463198935006588, 6.304815256362666, 7.119932713610217, 6.299855841806387, 6.587804285349651, 4.631973823146777, 4.56672049419411, 4.680205280472258, 5.316993407314629, 2.758412219465613, 2.734013847199286, 2.6685548442451954, 3.1187801022362853, 1.2644123728973442, 0.792792940730153, 0.02307475008993788, 0.0]),
array([16.36088576485368, 17.047810228795992, 16.508397925881656, 18.412431962609954, 18.519391187508887, 17.703680047384587, 18.49469589459746, 16.677667718971488, 16.28081030512009, 17.933373209625692, 17.33848682761915, 17.885101338102196, 17.124752514505623, 16.040937503425717, 16.126461696160902, 18.430845340751198, 16.831072123587287, 16.215716315345652, 18.14639847539281, 17.286506994521645, 17.053392965667157, 18.756250319080163, 17.12680198390337, 18.862734506587454, 18.83753082336623, 16.27432966334663, 16.93312731296, 17.72162179219033, 15.698658853749267, 15.863743734997835, 15.60351088598647, 15.554199171769099, 15.009096346676893, 15.24882837479615, 15.521110311677633, 14.791531737568334, 15.862987893759156, 15.76476913075436, 15.555465944507434, 14.994539575112624]),
array([17.989111181998613, 18.405851060970228]),
array([17.883642869275416, 16.438975949693447, 18.12610451209376, 17.726808701904528, 16.507705961562003, 17.067457760608995, 17.971756360712327, 17.009431443743857, 16.756364410005062, 16.80486548662126, 17.96525420967666, 15.724771246868956, 15.883942874763955, 15.568594640471561]),
array([17.008398082735436, 17.218144832283844, 17.435023582031818, 17.877046044252044, 17.91617071743287, 16.97534817989927, 17.667839666761864, 16.74268811118103, 16.77356717803926, 16.572924079713093, 16.627576757977614, 16.956745970267782, 14.660406861872556, 14.886094406330624, 14.406872267040082, 14.544304677282902, 14.20346592491499, 13.944052699581567, 15.219772800743959, 15.86990269611372, 14.946004998599026, 14.54971833441128, 14.174075357947325, 14.493383910312609, 15.01681248140415, 13.995786482639849, 13.950126285448055, 13.839799541524958, 15.181125504578514, 13.240462343562434, 11.72595886158284, 13.692658322686276, 7.8244800993255526, 8.462135433512728, 10.088615102525951, 11.47594335649108, 7.264982705537272, 8.12665056036589, 8.917094051125643, 11.41838984529927, 8.50236305423626, 6.9594730483128, 6.773053436959775, 5.866972257911188, 5.947492321491462, 4.055652226243548, 4.612356342172381, 4.013651618845853, 5.323294477843401]),
array([16.274835613142862, 16.73615856876533, 14.137128359941183, 14.567143828859015, 15.590985003531525, 8.693751742610488, 10.28402064890467, 3.7033271221773525, 0.0]),
array([17.81118968079497, 18.009360565137413, 18.132678185960156, 17.797881289432382, 17.71614054302086]),
array([17.23480304203646]),
array([16.31915359319142, 16.450844062967462, 15.958791342855172]),
array([16.864543429359365, 16.87929750799067, 14.863157022334688, 15.687716913160024, 15.057880539249004, 14.33170479968699, 14.2048759322067]),
array([17.34750672453575]),
array([16.803884005791822, 17.029036146833086, 16.884026721577793, 17.336515244899726, 16.8080342467614, 17.281531034846733]),
array([14.20638727001446, 14.63977035667486, 12.68954870275282, 12.908219427334416, 9.109028406597503, 8.579516806221779, 8.356464807920833, 11.313993267491867, 10.502707091427673, 6.866861441743151, 5.950651667364281, 4.922771520485961, 4.941729332253845, 3.845122314896434, 3.4585890194933344, 3.0295654157825775, 1.636977743534981, 1.2128354145888145, 0.21840220634362184, 0.0]),
array([16.326167263206205, 16.30808031588422, 16.107219505564824, 16.114817301256494, 16.294221580044788]),
array([16.378404919768528, 16.130442075570986, 16.25704272460938, 15.14137876811256, 15.299535601538983, 15.564050825295793, 15.259412447050329, 14.425142668677367, 14.879769766922776, 15.489804875696157, 14.03677480271715, 14.2402951699563, 14.289028509733466, 14.377574907906213, 14.915953365786677]),
array([15.998956028441489, 14.254253867644657, 14.772183770053468]),
array([15.975903667052249, 14.746393124694663, 15.727062143587773, 14.89930639289672, 15.42814066717917, 15.53131429369928, 15.574872625806647, 15.81479354229529, 15.005114460618415, 15.884152534555602, 14.48092253854988, 14.785452619045016]),
array([14.234162445995414, 13.92553522151087, 13.929990777418265, 15.378757408564358, 14.826281215911845, 13.974175119200543, 15.410460169882278, 15.9323132020097, 14.96659528612459, 15.358599388025262, 15.636279639562096, 13.911845713948528, 13.293575869477996, 13.585220829496768, 13.271925535428062, 13.675562887860472, 13.302608316778635, 9.630588016369249, 9.901732184977714, 8.582910964317941, 9.00854717675978, 10.696079848382688, 7.342620539706042, 11.238288137570054, 9.33576637730532, 10.4248799567955, 11.406732663790217, 7.684945123777908, 10.617179704960574, 6.347473788071432, 6.409601579370761, 5.984009464793113, 6.398567090836211, 6.476353695331597, 7.0603118210899956, 6.033825866175528, 6.682638028867494, 5.980243813113169, 4.20878722830683, 4.1911664229949395, 4.367304501421357, 4.370447270514594, 5.054563546712379, 3.307990584394525, 3.1563741962802276, 3.504438150015293, 3.1748014021163353, 3.3342048048873565, 3.034093045827675, 3.4468455653832732, 1.3934574216780304, 0.0]),
array([15.241412392114519, 14.151111923638865, 14.111662836307266, 14.178568342550818, 15.249589296861654, 14.362723163803988, 9.481812217858634, 11.02916635978658, 7.983708461215979, 6.2070736941923075, 7.156872051094779, 6.550817720223147, 5.936348750613927, 3.877779250282397, 3.7928715830374897, 3.536926468451574, 3.2247420923577943, 2.7579312970528345, 0.29961423998025133, 0.0]),
array([14.672198738584141, 14.724329921767564]),
array([14.567281969591319, 11.781953680001795, 12.167622968880464, 8.78399407639019, 10.823703030630739, 6.364794869451805, 5.373743801413773, 4.492279320083417, 0.0]),
array([14.337571282151456, 14.777834720714957, 14.641308090844133, 12.288803583675755]),
array([13.850111518329264, 13.829899634598867, 14.2992142186203, 15.053487091254846, 10.065692398177726, 11.49510981756753, 7.137520898909852, 5.994602582266859, 4.503624378784769, 0.0]),
array([14.443479419234855, 12.914442043302058]),
array([14.582752442219878, 13.37087274966956, 13.220414774880643, 12.850229428144525, 12.367933565070036, 11.89191909446387, 10.346077456369773, 11.3627785362626, 10.766721970840905, 11.488729627415674]),
array([14.372869601464933]),
array([11.734773336763283, 9.428307560985434, 10.223311776370638]),
array([12.795482267273512]),
array([12.206557179043346, 7.5365074366692095, 8.533076189453002, 9.812811922566029]),
array([12.317864057138822, 12.180265676860968, 12.07009688396257]),
array([11.657379142715524, 12.041142009667182, 10.960806788571434, 10.464168199207046, 11.52578642363367, 10.682798270459196]),
array([9.882724194061053, 10.950634281455672, 11.034407369125923, 10.824547893891323]),
array([7.904507465262349, 8.199461989449077, 8.573699697172483, 6.327695611822941, 2.5998240511471393]),
array([10.577394648889753]),
array([8.697121496006059, 10.474092778760756, 10.380439599694322, 8.490459906207887, 7.865991518736335, 9.949572307716334, 9.71603513454073, 9.292127866346478, 9.039621545259733, 9.690068694569845, 8.893809178130127, 8.140365359137686, 7.205306093483203, 6.942898659392661, 6.233233779094135, 5.625038302773765, 6.97262977258094, 6.638865483202931, 5.630967027326894, 4.717695744035743, 4.758435427260641, 3.8548641392865255, 4.607535784022208, 5.261469956894591, 3.556804555504203, 3.01763751172627, 3.180500121912585, 3.498759514802133, 1.9510346087614376, 1.3702230194958025, 1.3259888623389497, 1.3530139423572756, 0.13198339965690942, 0.09577555761492573, 0.0]),
array([9.313182702200129, 9.486858311726992, 9.33868814067806]),
array([10.318991697955706]),
array([8.86914277798087, 9.567096654634486, 8.628684602205187, 6.451638192938746, 6.1470890770173945, 6.471186318393755, 4.944109259104371, 4.658363703345497, 3.2123841494030456, 3.230087889789621, 1.7972091356708293]),
array([9.137915017764579, 9.334792339224585, 9.357876189545252, 10.059113443814391, 9.543083323976589, 9.320904042724486]),
array([9.861699491386025, 8.159564087814182, 7.3263798168921355, 8.715300242848286, 8.371601287819677, 6.637797299535837, 5.627019448721466, 6.365877509601418, 5.657066607968938, 5.685311419015835, 6.913427745506431, 5.4021575496630545, 3.8191500782980095, 4.993025076847791, 3.749475590177276, 3.7751229564377518, 3.0944675847758343, 3.4364141484643955, 3.309903512444149, 2.8264099657915174]),
array([8.650548844069736, 9.081887828668936, 7.256200949786582, 8.753561869798435, 8.699987331273935, 8.621752764052967, 8.56429997136455, 7.01764674714202, 6.235797120978618, 6.137357807235576, 5.973132422306646, 6.438755298703168, 7.121977522414692, 5.491742896272051, 6.319540814344648, 5.79899149841501, 5.067193406613659, 3.6024599573559124, 3.6462403854490915, 5.252348300993773, 3.690080607341166, 3.6998856114336407, 3.2127339797245655, 3.0339523943880473, 3.055885594265569, 2.761960194248947, 2.174032762793483, 1.6658349781057207, 0.0]),
array([8.448748904101961, 8.039589963768975, 7.66974838196491, 8.773796690715937, 9.680833340887823, 8.47654825473078, 8.247821859730822, 5.805803396823266, 5.536605512346174, 6.885785039362977, 5.488652115107057, 6.8971615598960465, 4.268456857472106]),
array([8.50409674307376, 7.448467004268305, 7.559869997471599, 9.87102020752428, 7.4309671399548645, 9.220026669315187, 9.307699534115633, 9.884794517817031, 8.19825408082941, 8.205463235364414, 5.52159062426063, 6.585663528047588, 7.08201508431908, 5.916848910254043, 5.4727992974926245, 6.836557789440889, 6.946590456382441, 4.3102707221783145, 5.195559751477307, 4.112087797859625, 4.344815910163658, 4.63239450017827, 3.8692056015766854, 2.66855683203499, 3.3161222819343528, 3.514801119197329, 2.982613846666811, 3.536547642247556, 3.425502247664487, 2.8453491049467545, 2.4415844574702152, 1.5065118294754205, 1.3764737615750025, 0.03152675052026059, 0.0]),
array([9.3411945225191]),
array([7.831478453042674, 7.283986356424615, 8.67661192978356, 8.146072213969154, 9.237799809444086, 7.7635244113231785, 8.776475174314877, 5.787384142457668, 5.606181502848524, 6.619581809101873, 6.1818871455078614, 6.449141220171357, 6.034036054000244]),
array([8.4420708165357]),
array([8.000898349053772, 7.533832603255876, 7.372843638002657, 7.7914048911564695, 7.305823464482174, 6.100037620234174, 6.620790604651488, 6.099937573049813, 6.6731577569794025, 6.548887671833928, 7.144353466956439]),
array([8.071479703590626, 8.519401379839532, 7.471874501138046, 5.981165287160192, 6.075273777187862, 4.836055389387401]),
array([8.497621362346727, 8.006966121889016, 7.753018217745652, 7.05834420411627]),
array([7.429106933830127, 7.841815347078685, 7.728467621933358, 7.798362977656569, 6.69812730237631, 6.073511923487115, 7.22648380924094, 7.069314960326816, 6.425837848117513, 6.466315364299615, 6.078505955847259, 6.939563317645537, 7.175605678206019]),
array([7.230422722193796]),
array([6.954058459673786, 6.033933921023284, 3.6625309432421895, 3.0935295557735616, 2.0222383000871456]),
array([5.912072027939607, 6.863270599281237, 6.0661132605419255, 4.397837472982504, 4.1241182223659525, 3.8309360929619807, 5.12633233578562, 5.116025371447403, 2.9314864793153936, 3.392469024481441, 3.0098375779563904, 2.782779139696907, 2.645523449841772, 1.780386095630134]),
array([7.336732272403632, 7.227647397808471, 6.393804003949693, 6.625898806863654, 5.496674410592714, 2.7615493149043004, 2.5455555147124467]),
array([7.466629042872865, 5.799130464590062, 7.064204383981401, 7.021833275298737, 7.157947388221173, 5.755113703830117, 6.433941089315606, 7.112988823803806, 7.105448046919658, 5.605742479921989, 3.6892836866723213, 4.844808390995287, 4.8076937228106145, 4.260106373759994, 3.814800075405463, 3.973471007891235, 4.512224830714609, 2.7663305529256332, 2.217908673897539, 0.4272089413961589, 0.0]),
array([5.90511932949977, 5.568376936107487, 5.866263614508373, 5.766121163501202, 6.803992752761693, 5.907204379944529, 6.8191995244598465, 4.992483750452123, 5.045695426265346, 5.304346882497324, 5.283321067596125, 4.862794585835711, 5.205462704072234]),
array([7.041553741358552, 7.112637393097513, 6.588243172789969]),
array([4.178449686121207, 4.1093160281985, 2.858457666008271, 3.2090061931628604, 1.262430636202411, 0.0]),
array([6.359680530510662, 5.850235617918439]),
array([6.0012294449310595, 3.931834801779777, 3.3741710905096416]),
array([6.01823505221126, 5.7982158619349375, 5.963620720158826, 5.81461162586353, 5.690208989834678, 4.786140461538421, 4.783825343939925, 4.642207346411053, 4.224932659456125, 4.151872615346065]),
array([5.846748028688912, 5.644485578994152, 5.249694965041585, 5.052859994021661, 4.898521232860337]),
array([4.537041559442616, 0.0]),
array([5.751186307281383, 6.034080261876153, 5.095000675096936, 4.736961233020776]),
array([5.360184714213317, 5.962377448279686]),
array([5.746018203119284, 5.792066907573569, 5.427687982510959, 5.515683241134868, 4.979391975374686, 4.837944238444429, 3.109846488056026, 3.0674338864241304, 3.4037194751015183, 3.136771428621019, 2.5914437508969264, 2.4414201630756915]),
array([5.519629576417434, 4.681385646347736, 4.681288622264352, 4.979149848778484, 4.612431897042526, 4.285792817130982, 4.786251166218825, 4.565915523362179]),
array([4.004462051085701, 3.756688268105629]),
array([4.546296995588271]),
array([3.69634266443867, 4.9467717788090715, 3.8723900354223977, 4.374904942423651, 2.896140619306066, 2.85477575271844, 3.2032072261861098, 0.8069722078274536, 0.04921429485023704, 0.0]),
array([4.42820650282634, 3.73791664332209, 4.483066317244642, 3.6406205029733165, 4.278556055987294, 3.8921902909689585, 4.730975031763281, 3.765375147508285, 3.881725588518995, 2.881309340953323, 2.9432987133739528, 3.1211687673913735, 3.509425657801273, 2.9967038624328755, 2.7659525437669226, 3.1121320496332596, 2.8303455113466334, 3.3870756170243896, 3.1090571238456186, 3.3574598014797767, 3.244369560860733, 3.4800685112905554, 1.8797248363108787, 1.5648173825801246, 0.8611221922470853, 0.5038384532192266, 0.5959111910932928, 0.08366784440452554, 0.0]),
array([4.764961991823351, 2.800941087748716, 3.5781369919135138, 3.0437468031693022, 2.283923451089, 0.48271137059598807, 0.08121400167511211, 0.0]),
array([3.8081173148315646, 4.470981938653432]),
array([3.6259454145881316, 3.844639067002885, 4.4378384930728645, 3.8166052452848125, 3.554113903335194, 3.1449449753862124, 3.4586706657777997, 2.8158810964306715, 3.4149530637753136, 3.06664598477935]),
array([3.695463631921913, 4.047662318454358, 3.0497700371068914, 3.2108525022192613, 3.483671198591785, 3.131905650698022]),
array([4.058029214510696, 3.619675008326714, 3.4057248664335464, 0.9833479249708577, 0.0]),
array([3.7668892340030853, 3.7998012142010507, 3.76401740835617, 3.6064043621844184, 3.7452552782276367, 3.8007389465825043, 4.238385223667438, 3.751883606965399, 4.226090132064765, 3.887867681550952, 3.4347460726651553, 3.1906281559655163, 2.7353682698248925, 2.8294886838834716, 2.9636041571660203, 3.387799563717464, 3.113714835082515, 2.8607315762679444, 3.233686031209593, 3.2473278643122843, 2.750219929221527, 2.93159827466275, 3.57271056688788, 3.547644498135836, 2.591885013923079, 3.2876860429980312, 2.741770369982543, 2.7426702853144915, 2.6997628269916945, 2.8032856847074417, 2.8736013229968442, 2.010487748203585, 1.8529976562086032, 1.6725338720768654, 1.641022813865734, 0.4283704032355812, 0.4614004959662695, 0.5931046160439492, 0.462572613773085, 0.4560334939455346, 0.0822789901268566, 0.0]),
array([4.246833866335378, 3.964061222862449, 3.21973406105743, 3.078233412384014, 3.049099813144391, 0.4412317714419655, 0.0]),
array([4.221365767810085, 2.7075462603399147, 3.329809327862656, 3.126755809136401]),
array([3.669184004625472, 3.751971632959558, 3.785548840762692, 3.536810789895473, 3.3445631616094813, 3.5396844439255992, 3.0090010309774566, 2.734603906790396]),
array([3.9096724342380758, 3.689366035114409]),
array([3.121158319363605, 3.5986224297227944, 3.3377655658585317]),
array([4.105042624739935, 3.6749978792029667, 4.051821545722468, 4.173368505324024, 3.0411387249129884, 2.9385795696117722, 3.295167566223199, 3.1579839768320337, 3.086720248052937, 3.0862190901898834, 1.4615282789269117, 0.9194434336727312, 1.652253807649771, 1.6114981808697881, 0.1879320746494232, 0.06387117405304253, 0.0]),
array([3.727262198607161, 3.7355473862151314, 3.6475675819672624, 2.9291751990299475, 2.9388106368560476, 3.235795725915705, 3.307687195844387, 2.948226743923953, 3.0209113908274645]),
array([3.6149217873894663, 3.908942409313534, 3.913615026145013, 3.9126178415380055, 3.2346465721258477, 2.81550170783819, 2.873100519626188, 2.2803289765792187, 1.1898454832385736, 1.6372983334738325, 0.19508555578319275, 0.0]),
array([3.7043266813678044, 2.901338995198845, 2.898276269129526, 3.557648786373993]),
array([3.0687499833055436]),
array([3.044566657251629, 3.1255758902110524, 1.3888813406556757, 1.0946357675963494, 1.327525501907417, 1.172435233710981, 0.5394457675531709, 0.01521051978157223, 0.02650168107341551, 0.0]),
array([2.9519888888684336]),
array([2.785930193115341, 1.7511356969848886, 0.11632572651944197, 0.0]),
array([3.058338792471935, 2.8780904105767386]),
array([1.8902627820284996, 1.4800596774226973, 0.9699636903167274, 1.4241971778612006, 0.9402670877217811, 0.5641559103196729, 0.2214594762251818, 0.06039907590244878, 0.03868140664250544, 0.11017423426491693, 0.0]),
array([1.9621304752452713, 1.4499717458775052, 1.284052339510201, 0.09494837072745602, 0.0]),
array([1.350087893863981, 0.5023875160613029, 0.03555059707857651, 0.0]),
array([0.8016727426276757]),
array([1.56797516032673, 1.0703089868366091, 0.0]),
array([1.9800021333176099, 1.1353911648150528, 0.7027941559122576]),
array([0.07468355038812115, 0.05030956515843324, 0.0]),
array([1.8803335652135378, 1.6858012338953272, 1.53825091939791]),
array([1.1254550509227212, 1.0557550686632093]),
array([0.002381059000615096, 0.11114755236578394, 0.04873782404682049, 0.0]),
array([0.01663370248756592, 0.0]),
array([0.39397028864322087]),
array([0.7319662465754379, 0.0]),
array([1.0883585072087056, 0.07588124594406728, 0.0]),
array([0.39031922519655127]),
array([1.049222851827658, 1.2662546227924611, 0.9487157537680966, 0.499404011279234, 0.1255155841693265, 0.0]),
array([0.11503458003698114, 0.0]),
array([0.8013241062012624]),
array([0.47173799191274174, 0.0]),
array([0.0014782385461840075, 0.0]),
array([0.024175155812469637, 0.0]),
array([0.050284552888078135, 0.0]),
array([0.02772697478199236, 0.039073654369542674, 0.04940871688607526, 0.1151524389027829, 0.0103706663651673, 0.0]),
array([0.10435250543364782, 0.0])
]
d = [data_1]
names = ["13"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T8', 'T9', 'T10', 'T11', 'T13', 'T15', 'T16', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T29', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T42', 'T44', 'T45', 'T46', 'T47', 'T48', 'T49', 'T50', 'T53', 'T55', 'T56', 'T57', 'T58', 'T60', 'T61', 'T62', 'T63', 'T65', 'T66', 'T67', 'T68', 'T69', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T78', 'T79', 'T81', 'T82', 'T83', 'T84', 'T85', 'T86', 'T87', 'T88', 'T89', 'T91', 'T92', 'T93', 'T94', 'T95', 'T98', 'T99', 'T100', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T108', 'T109', 'T110', 'T112', 'T113', 'T114', 'T116', 'T118', 'T119', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T129', 'T130', 'T131', 'T133', 'T134', 'T137', 'T138', 'T139', 'T140', 'T141', 'T142', 'T143', 'T145', 'T146', 'T147', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T156', 'T157', 'T158', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T169', 'T170', 'T171', 'T175', 'T178', 'T181', 'T182', 'T183', 'T185', 'T186', 'T187', 'T188', 'T189', 'T192', 'T193', 'T194', 'T196', 'T197', 'T198', 'T199', 'T201', 'T203', 'T206', 'T207', 'T209', 'T210', 'T212']
def get_taxa_names(): return taxa_names