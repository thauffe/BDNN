#!/usr/bin/env python
from numpy import *
data_1 = [
array([29.254213358271542, 34.12746729264393, 30.49056798426628, 33.8393451474678, 32.98337156578529, 30.659071677782336, 30.447510230180356, 29.53240432890513, 31.865513863280665]),
array([30.911442834720404, 30.530029064008254, 30.473783492686806, 33.30164745707374, 34.26259391267158, 31.977223981648912]),
array([33.40552028994163, 34.38375925109043, 34.022580430150484, 32.93147172849999, 34.383110678327085, 33.601283873139224, 34.297198895363195, 33.115907333254505]),
array([28.811599412551097, 32.988470656648715, 33.226175429056404, 34.07937901364313, 29.551613918198314, 31.480431413552584, 32.02683828236445, 28.20362981134082, 33.82345606448176, 31.54841787733731, 32.86348846230157, 28.80225577786573, 30.99988620018852, 31.270413679180244, 32.172995428624745, 31.33285385143751, 30.290072649583337, 32.82497893100037, 32.71439877933689, 31.33380646670123, 28.847061889981283, 30.16486761530508, 30.9483322201969, 29.109538729558132, 32.52162231810503, 33.175712203123666, 29.100702369435467, 32.168063502999345, 33.612255466124026, 24.00660108876263, 27.969656229428995, 27.623238105675803, 25.226306328412306, 26.57379491245648, 26.816795231250083, 20.529399244035424]),
array([33.63959740740995, 33.35187878046638, 33.556833718219, 33.156135514392915, 33.73631917809499, 33.43909403697191]),
array([33.53764403550919]),
array([30.383966659961043, 31.147706837379225, 30.138014782369467, 30.411881947499662, 30.9935316415968, 30.814384053501726, 30.940887651120946]),
array([30.761404873876266, 31.541389798865428, 31.501639383786724, 31.13959306967421, 31.06005872743007, 31.869054219496988, 30.91978579196099, 31.44489465997008, 31.77616380180302]),
array([29.32590179309965, 30.371580391987933, 29.5560279633711, 28.66621138400827, 31.581909551315235, 29.987710421991615, 30.61402186151375, 30.128444262975826, 28.29125854281019, 30.657012617872475, 30.1146540647516, 28.367309868022197, 27.42592004858797, 23.69361320707958, 26.05540750951606]),
array([29.620980564939025, 30.26318368903994, 31.082360638080992, 31.28412136314645, 30.053077505952523, 30.45471143190295, 30.856522630914018, 31.21583244994114, 30.702901366577414, 30.549034787042462, 30.84396871665585, 30.957589446519652, 30.22607805992815, 30.982986822415963, 29.89426297575365, 29.932134956599675, 31.2064397047646, 30.05475949617416, 31.05905878301636, 30.666966629254105, 29.816136221234324, 29.73266772707253, 30.450395493671728]),
array([30.49164597804767, 30.543687983543553, 30.3102893168199, 30.312303944199943, 30.384325378819078]),
array([29.341513577476817, 29.61085974366924, 30.313728195033416, 28.71466799306188, 28.19491098457553, 27.08315252376822, 27.212817702616373, 26.24518280286035, 20.711588144250005, 19.943643293223715, 19.44152142528189, 20.000459647206256, 19.06216374838123]),
array([29.159018076722642]),
array([28.95868964657858]),
array([28.336889268309207, 28.840044524046437, 28.739216415677404, 28.33669100771033, 27.531407356918116, 26.861668221528443, 27.94061815361103, 26.973186364531887, 27.515462348670603]),
array([27.812837494271996]),
array([26.170083193734282, 23.19730036109117, 25.20793579735049, 25.566783256630664, 26.248905159257184, 24.239142722037194]),
array([25.834215427501267, 22.721215095173466, 22.173916093955683]),
array([26.716940195782147, 25.71979152980959, 25.558839835381733]),
array([23.544056738612248, 23.888555965590296, 24.143864683491188, 24.05845351883514, 23.235702263843706, 25.82990994865427, 25.10459691828631, 20.600917778421707, 22.584732890332603, 21.311141109136585, 22.046564205118724, 20.50521818120645, 22.611414206680745, 22.075844856268844, 21.39190830200556, 18.428085552783138, 19.08775057674107, 17.71478962685759, 19.77025304843151, 17.684573862970076, 17.63034402906376, 19.913635552849655, 20.38787546291208, 19.38870214694694, 17.953821888559162, 17.384549130562437, 18.957717779187583, 20.082330004964227, 16.760164566415767, 17.66529586516704, 17.090562976611675, 18.561952821965715, 20.110424794451166, 18.087618209112357, 16.7632085142042, 18.870713861029028, 19.597665563055674, 19.4322001678594, 16.83406344914454, 16.76644587140114, 17.835058132760516, 15.975766677717537, 15.994784082558061, 16.852969462643006, 17.33728949328385, 20.351817899965816, 14.765154549270807, 14.905195971200357, 14.228527497979385, 14.191154230813432, 13.98852995451093, 14.196880151801741, 14.291094228850572, 13.04221484120184, 13.194784288317141, 13.082519357710021, 12.844930465545982, 13.126521039940087, 13.380096081712951, 13.328600301576989, 13.196836554338125]),
array([24.36940226374295, 18.82838585847467, 18.0766241096102, 20.2443746119069, 20.246710665383603, 19.837361538803187, 19.406455080448804, 18.695923387483838, 19.875292739464378, 18.941904533013997]),
array([23.99289075727824]),
array([24.831438375835727, 24.64905226794904, 24.23490626000369, 24.367826557176613, 20.394040235695222, 19.166890035934955, 19.20959384387594]),
array([24.115863858892332, 23.251503979624722, 23.655048258201695, 24.640079304604726, 18.90677833951647, 19.902224701681078, 18.405153451324683, 17.326156044386487, 19.876601095753962, 17.417056068045355, 16.679016476579026, 17.287214883228238, 18.12498393403937, 18.774780367578444, 20.213318214634675, 16.979767897819745, 18.360924158462982, 18.671273308480803, 16.16948168808054, 16.865730242517753, 19.100716597577787, 16.284458445345255, 19.386097135973156]),
array([23.440020013108505, 23.86206321640696, 23.130174883477473, 22.28624577282386, 22.964325887414955]),
array([24.161583092586547]),
array([23.535840870670903, 23.103497644615974, 22.180185986371026, 16.227059366843037, 17.170416608798316, 20.400401971301488, 17.774881354602275, 17.677671834891413, 19.189707178972572, 16.72957866949003, 17.507344194604492, 17.005579268710903, 19.68883373712506, 17.044106874559702, 17.802204031573957, 18.578491134547157, 20.1684982891461, 15.968422746697176]),
array([23.07345836792644]),
array([23.183398759126987, 21.76852747305714, 22.38004609072449]),
array([23.199346878292342, 20.547254544865666, 20.252539565924817, 19.636795690380893, 19.37388438542286, 19.5400447256525, 20.233282336201007, 19.34391674936593, 19.50326689388989, 20.00330564171458, 20.425938490243304, 19.11931824099708]),
array([21.10494091061072, 18.819551089144916, 18.4124727759755, 17.999565195827838, 20.11865625893331, 19.698243128185563, 20.081363140167724]),
array([22.613254842772633]),
array([19.53653980148772, 20.035748433431262, 19.07535377266873, 19.643680254975596]),
array([20.608568481039164, 20.4359643200518, 20.030316361304816, 19.974032502732193, 19.91487747759275]),
array([21.745989599196854, 21.758593364402973, 21.785713686158065, 22.15028640068873, 21.884886145694704]),
array([20.633649806378955, 18.625686040930095, 18.717102274931012, 20.231607163160767, 18.74530957426685, 13.407007629997148, 8.971388827887363, 10.18045208041476]),
array([21.516671465737726, 20.413966515944686, 18.137652144175036, 16.937179452744353, 16.272267987590684, 18.789574261610245, 16.088846664583798, 17.73510428762914, 16.507368823922036, 16.00339421594608, 19.682095922169587, 16.191750351440476, 17.79253225444781, 17.472464942902953, 14.7183723015792, 15.043700042997694, 14.796348909517075, 12.207928616462866, 13.23599654035248, 12.394902436204555, 8.77419754401377, 7.334280677445723, 7.560341485467131, 8.061966041745155, 8.640044035666172, 5.7109980088884855, 5.45985677895167, 5.729042133346693, 4.397810055080553, 4.752516344980978, 4.066279227604891, 2.8430264212999887, 3.014906482371016, 3.315504701321173, 3.501512488692098, 3.4136611229318654, 1.8490145424813744, 1.6077204655149668, 0.0]),
array([16.996101259363726, 16.738633563717524, 17.49620139749967, 17.00430262871405, 18.296774629178955, 19.68046841020802, 13.997993368935411, 15.220045494146019]),
array([19.294409079481305, 18.581843511839246, 18.907282508329914, 19.292896567991466, 18.898266486736514, 18.575125403727473, 18.63114848563782, 18.037967180517775, 20.051782643201992, 17.884556664912992, 20.02903382593258, 18.194618288276075, 17.74294842978679, 17.59710061323459, 18.581728203439827, 18.487587715870998]),
array([19.939609199267498, 20.219229521135937, 19.84725343755975, 20.01388911659122, 19.140753085218535, 19.211008428500257, 18.822808176615524, 19.297043732488444, 19.532622300351424, 19.249290028355187, 19.201267657109366, 19.070094705826683, 18.96257973954483]),
array([19.507916421021513, 19.586587901825897, 19.8817034519306]),
array([18.04351159079051, 18.220643905341294, 15.600583691652416]),
array([17.675769514871295, 17.19097219118973, 17.709214126952524, 17.0681262674916, 17.354252166017226, 18.97602968449241, 19.353798820116623, 18.19901516162238, 16.80188262509591, 17.519751408013953, 16.49612687822331, 19.525850770591738, 16.58387928379784, 16.227504888161768, 19.035425159863912, 17.06927252018237, 19.084505877963483, 19.38427801552187, 17.67044691509948, 16.42247655246609, 17.652412742206742, 19.321628891574036]),
array([18.550097024180275, 17.98969184430996, 19.851808092775585, 19.786914904130946, 19.129584436982828, 16.703877931887558, 19.673415517260224, 17.008367406710306, 19.274541167098686, 18.40041722848691, 17.786354253344992, 15.50918612558663]),
array([17.257640983634943, 17.56356041681864, 18.8624624724647, 17.981207292661026, 19.54812846784774, 14.643154246706262, 15.931863776544102, 8.466619513362504, 7.836153651007721, 8.462607620227494, 7.409922363471963, 7.322761305382364, 9.485254818467412, 8.926613457589035, 10.855115608774275]),
array([17.81998960103777, 18.27248682968796, 19.391576416381916, 18.02764796585811]),
array([17.321930258836556, 16.541501484284932, 18.44294264110855, 16.792096590583455, 17.231526092827306, 19.078959084149783, 18.791904758171484, 17.340059096994498, 19.230662515012565, 17.33276722262921, 18.005372900133608, 18.559682012403986, 17.30161623973755, 18.99397768091035, 19.202925386576347, 19.281282045348686, 16.386624022327716]),
array([16.810295470241034, 17.871994271420103, 16.04050416863533, 18.87396785939183, 17.314848814611267, 18.530883485944585, 15.99239785506519, 16.60558678338484, 18.482510759406107, 16.058164907789664]),
array([17.645897261743578, 16.9446819332602, 16.974317986856537, 17.238624241372726, 17.749598611137042, 16.77413072573323, 18.105810971950817]),
array([18.27716278583968, 16.797475553127246, 17.15161825931574, 16.15125040363126, 16.251163399351935, 16.47090233771513, 16.858358471177247, 16.304360098227786, 18.228440319164708, 17.081696994000392, 14.086382280435604, 15.329593953834847, 13.47296826852772, 12.227156519015999, 11.68632996214791, 11.033830773046274, 10.770648350962913, 11.024774152794743]),
array([18.20403172759045, 16.04535699047633, 16.990221122460866, 17.702386693863698, 16.773357081994362, 16.75309210858946, 16.058878197612778, 16.281538748838027, 16.22252847719426]),
array([15.830432369564281, 14.685946947040259]),
array([17.819389529665496, 18.11803826119353, 18.044585117938897, 17.645854620985386]),
array([15.714400189477985]),
array([16.47107974296307, 16.379708072834944, 16.121511546799738, 17.219495547604726, 15.685020194571486, 15.872395381026825]),
array([15.981427927297908, 16.219734473507756, 12.893236471306441]),
array([16.95522989197989, 15.971460827157888, 16.121490808052812, 15.457896738356878, 13.932500161368797, 12.708062281794279, 11.641091210556983, 11.794231960272109, 11.837086916231277, 13.312441080858472, 12.243555774866344, 7.316567951931938, 8.126508610442482, 7.8286429953990435, 11.146014543637591, 7.376915262755087, 8.408837861290811, 10.572251141988387, 7.59421213864095, 9.854903292247098, 8.037462846217608, 9.075022936159995, 9.436792282498056, 9.771103197370826, 7.747979573479899, 10.564086990095031, 9.430913581304786, 10.413552270786948, 8.835884972102326, 6.6804507846038526, 7.173786048291954, 5.938360079374761, 4.780928586221517, 4.437684884380905, 5.319686682836271, 3.55870023645418, 3.5140617140497765]),
array([16.7448962148613, 16.38065639876002, 15.276361585296922, 15.15819666853825]),
array([16.529845153838345, 16.62166734803582, 15.363701450123225, 13.231008579087659, 7.433261218949372, 7.024464004022882, 7.0857037358877735, 5.330300529350934, 6.241234554695563, 7.095670491112625, 6.658663406617192, 3.80140679878786, 3.8534767048744873, 5.029466514297694]),
array([16.086656064004433, 14.118477605449073, 14.7793104148351, 13.42258486412429, 11.969826413543416, 13.00263785571103, 11.886391517413738]),
array([16.243098002094545, 16.16760413752079, 15.779555228530027]),
array([12.784487495569945, 13.370583682287556, 10.71519004491165, 7.168604766834552, 6.58734766088642, 5.049827102231983, 5.028966186768485, 2.7963452864768676, 0.8928945963801123, 1.6921445307379819, 0.0]),
array([16.035223167488564, 15.007793249766491, 14.155919640649184, 15.409944269408687, 15.542915897390447, 11.945126014289166, 12.058469173270666, 13.179966759325392, 13.02666237812524, 13.3726479214096, 12.13065404736055, 13.528509085039538, 12.928159569698412, 12.301226823778283, 12.925986917341806, 8.29490341620679, 11.1372592727318, 8.569407160506685, 9.891076208270647, 9.971966547628707, 7.6821344470288935, 10.240554850460988, 7.548827371062491, 7.596570587721115, 10.481774010526411, 9.630645463541926, 8.443500518108328, 7.585500618842836, 11.50684610813202, 11.54380055688736, 10.430684012162846, 11.350526808166125, 10.495895322747222, 9.040138689051101, 9.74855114029557, 10.625345783974426, 11.194451238162808, 6.672945266646483, 6.869110586616024, 5.914448028404515, 6.263626265574502, 6.919693460821182, 5.70572623277837, 6.181244766398377, 6.559056416542306, 5.514322431583196, 6.0332803027299855, 5.506211693956758, 4.243404979076888, 3.631691409095443, 5.031574366420312, 4.368599953079077, 5.1683580749053135, 4.887230865358591, 5.151103789603782, 4.840729140613492, 4.797894638781529, 3.7592330777396046, 5.306610817949133, 2.6779194470732044, 3.1240599738583867, 3.0740174737378734, 3.295998916899432, 3.287999364617824, 1.9742316019100499, 1.0517760352991694, 1.7844088033501009, 0.5820451410402059, 0.0]),
array([15.463650820326682]),
array([15.330955930778583, 13.807287988698208]),
array([15.241530287918042, 15.378288327675413, 12.623824131121843, 13.280134007975182, 12.638235134117096, 11.779744785318364, 11.560243690567571, 9.799304491783596, 10.481151326229845, 10.280106267689344, 9.848084022230625, 10.984685354047574, 11.122147705140664, 10.92061288869182]),
array([13.286081965437855, 12.860946603876908, 13.350331195763312, 13.748456783494177, 11.007161725913203, 11.353475240382291, 11.484915868686032, 10.963856901232504, 11.28872768504818, 10.299950563438545, 9.893840131585687]),
array([14.638815981732009, 14.477225231583091, 14.284321113066497, 14.125610375170346, 14.64299611714417, 14.525635263622291, 14.987753543237638, 15.27399814905975, 12.9976349773331, 13.556620616547203, 12.068413901205387, 12.461521196004197, 13.326177825674725, 12.13164942527611, 11.721602466240821, 12.960462041789523, 12.212496696316979, 11.671895078855478, 11.720057023426138, 13.574122748494544, 13.1645894424325, 11.127162820125953, 10.985636062607002, 10.902659333985019, 10.762035304580005, 11.295177361123837, 11.134127780086278, 11.003789973828006]),
array([14.565451312318029, 14.834631263132826]),
array([14.444130864715815]),
array([14.127214464383176, 13.526747721152999, 12.378615388808011, 11.46950782497257, 7.796888725737528, 8.095242090547076, 10.319126013264608, 11.36751225992653, 11.242183387546866, 10.394127980298562, 11.082253970886004, 7.848420960689809, 3.659338273881059, 4.806702087183251, 3.333686287429987]),
array([14.947162590400936, 14.77308782872153]),
array([14.668254023284291, 12.727272290430964, 12.988981368263666, 13.698847167093112, 13.723595710583668, 12.795199138675866, 11.596686778402551, 11.052754323043878]),
array([14.21040135043238, 14.573076477721534]),
array([12.330660016298218, 13.131229221464363, 11.589128030966076]),
array([14.214233646499734, 12.998556125286774, 13.803478622217238, 12.237963549453063, 13.323651228785527, 12.570420180663296, 8.086311106095605, 7.9462738955139764, 7.900323224397928, 9.150009942466074, 7.211544156688114, 6.141805639300686]),
array([13.915057021805573, 13.563207666891158, 12.69696816856181, 13.342790579496567, 12.999611573105582, 13.755205770881597, 10.386286784870196, 10.900395139425978, 8.89051506511776, 11.414609993604277, 10.225124317113321, 6.520738342902865, 6.836571284078533, 5.791713507702387, 6.033398013391042, 6.845009232297768, 6.554515480427605, 6.90196618634797, 3.810974338037232, 4.977920752992123, 5.1456693076342095, 4.231483215051236, 4.697945536566009, 3.1320661717266516, 3.489194127162103, 3.356127950437253, 2.412930683853394, 1.7928215151969145, 1.164257378416643, 0.5193046464342848, 0.23127661313821657, 0.0]),
array([13.99804311012577, 13.601406518184394]),
array([13.181079994699148, 11.118725037738674, 11.624982719657739]),
array([11.99805525745796, 11.7777633063978, 13.4123999657262, 13.456475021469476, 11.2757733573535, 10.46046699596903, 10.430270069626053, 11.413547957145717]),
array([12.594672871116183, 9.131164407985132, 7.495789600722234, 5.734858872936456, 7.1778290316825615]),
array([11.284860909971869, 10.334654315995607, 9.334413253313254]),
array([12.06205063500766, 11.549103438736843]),
array([11.900918192169229, 11.83556624765215, 10.883103557173722]),
array([11.881669035516975, 12.155548636472643, 10.466111328174692, 7.357035075591772, 11.208519013801668, 7.406858313877866, 10.572853077015237, 9.852639954471986, 10.691498567514952, 8.916606940825625, 10.606105813550124, 10.665366343998397, 8.840161641932639, 10.21237187410712, 8.099343100927742, 8.443845920975667, 8.340658494810144, 7.821885730832706, 7.377232132996179, 10.287127752206635, 7.971610157066548, 8.05515251527313, 9.447809496245572, 11.092268113983266, 10.806533660184753, 6.6735520427831085, 6.5603692415844055, 5.845905559109157, 6.7318323458443805, 6.545837662370138, 6.7732708358860725, 6.056687032894509, 5.387792198224492, 6.647704254004444, 6.621956244272248, 6.958472764589267, 6.739050009543519, 6.13911222025691, 6.697395439379322, 6.713356917907529, 5.855389431830777, 4.616182349484675, 4.597931516071242, 4.7391958452471385, 3.1950904077347513, 3.3930827970905346, 2.6343673650142403, 3.5486959823172186, 3.3480386945345177, 3.2548841828417614, 1.9713141310349367, 1.7645750627695782, 0.8248997769616337, 1.5304133696595663, 0.9937493904909704, 0.8596277657959958, 1.7619726060994456, 0.9909911177497126, 1.341756226195889, 0.6720159268550279, 0.0]),
array([12.401710701640877, 12.341719476923513]),
array([11.74913659523872, 12.023996227701216, 10.13087224191654, 7.265098705140557, 7.321147942819825, 9.2652864218698, 8.028091859738135, 7.796330832818526, 10.371955101884797, 8.058111022503955, 9.753951227315225, 9.004012731132203, 11.100195605250958, 7.366183134886869, 8.623033258262394, 10.412713406419979, 7.960957979938213, 10.228637035556309, 10.591616863955156, 8.720897584595843, 7.6455546899094475, 11.355873039229278, 9.354635939932267, 11.226600592714018, 10.207122897877966, 10.313442235670918, 9.328815866969556, 10.369576750632703, 6.823057824639604, 6.978504460397345, 6.53580736251506, 6.899927500889781, 7.146949201081476, 7.232842372542125, 7.071207291576172, 6.286576851711667, 6.1715784571796934]),
array([10.681434012477512, 9.752769540795946, 8.201343225328358, 9.007362617662315, 7.841841611663174, 8.076534709093906, 10.41504088283759, 8.260031324364618, 5.785188488058392, 6.969973501117529, 6.3060169746885055, 6.468068647274668, 5.4799016231226005, 5.610108819448438, 4.722305113988887, 4.458373732171248, 4.083258277805958, 5.199660859290772, 2.831893153516345]),
array([9.736370423908012, 9.820316076845987, 11.341242334574556, 11.573821560143507, 7.74338308677909, 9.709698471496996, 8.357040776743915, 7.534775996386665, 9.818356625385986, 8.707552562687997, 11.246485342058376, 10.514515590007235, 7.502549700076118, 9.98033116442404]),
array([6.784615063824765]),
array([9.908272227035912, 9.13362941364845, 10.21235732313584, 10.814745969050836, 7.804249905626069, 4.173832724565015, 3.361036078820159, 3.276812596857946, 0.9408630299858461, 1.4053124401522379, 1.3297357889654398, 0.0]),
array([9.056613729380224, 4.003723491909091, 4.835755825911553]),
array([8.41941304414844, 8.746310567778561, 10.197342240567941, 9.205252489154118, 10.068330406605558, 5.819075078435118, 7.093403752261444, 6.441883763483269, 5.507319372249301, 4.208847859283281, 2.9820272875920626, 3.4834048942920015, 0.0]),
array([7.859775337226891, 10.446483626267012, 8.811712406132855, 10.931210836856756, 9.144528684507268, 7.827173206580984, 11.196735253162972, 9.928347267369995, 10.552835522949408, 9.83508474557891, 10.675490978641914, 10.438718333378743, 10.902233333527668, 10.067279532199542, 10.893858308590227, 11.098741602423454, 8.703401606756236, 7.529724283519535, 9.760011461294434, 10.728486706079995, 5.476849564406255, 6.845586048341033, 6.311387894464991, 5.373307821884938, 6.101173872899326, 6.334037056901054, 6.65733397615447, 5.814071531726128, 6.0569417733029844, 5.702966953919317, 6.369115159427336, 7.10849167442609, 7.039904612594031, 7.073507721708805, 7.129122144164639, 6.7377082813682065, 6.545202256652146, 5.285527156403561]),
array([10.821363118670307, 9.074284454980035, 9.593611406600695, 10.39664606006459, 10.657282227134258, 9.17435565603081, 9.779144644138071, 9.971119364741861, 10.719447796498702, 9.083299715306879, 10.53083466808371, 10.172507434287251]),
array([9.222689088393121, 9.048357185232739, 9.363775605673366, 10.827294366946886, 8.91976470086937, 10.476392920358997, 9.471787848083352, 10.162103668352676]),
array([10.623567687266013, 9.869351390126548]),
array([8.175755840479177, 8.216219223035276, 6.3335715918082665, 7.038821026902268]),
array([9.366069503461723, 7.310346550272227, 7.59393151473791, 7.459232455671538, 6.14411064865061, 5.555526533561111, 6.123383992944813, 5.455038769496937, 6.857434709974319, 6.010431951348499, 5.943768293672207, 5.642164363508962, 5.3558212717520695, 6.669263322988763, 4.031780087704933, 3.8675625803777356, 4.747891891956725, 4.377699244685243, 5.038811563289136, 4.1213328302798455, 3.603306561430698, 2.7540896961597205, 2.1906296627123436, 2.276629287789725]),
array([9.334337297096875, 8.715401044840979, 7.51377854016183, 8.999593341686586, 8.308144877909678, 7.849801042054104, 7.809502242510258, 6.823512631970205, 6.934908289426596, 7.0176510421032265, 7.208792441553186, 7.1470250718978265]),
array([9.071661720321584, 9.243896275370096, 8.131562171979851, 8.038485439150804, 8.235734796723607, 6.961866632756579, 6.156932614798997, 6.8040830749919445, 5.666147069088203, 6.798572651199904, 7.197866264968297, 5.9481774595474555, 5.65588033033341, 5.545430515009107, 5.968872486222252, 5.516208769177725, 6.8623016893346165, 3.6722538395686652, 3.619604256000991, 3.934572342072032, 3.7129696676005173, 4.599589268038298, 4.131656557143868, 5.109357049464167, 4.184735052313732, 3.841834796409649, 3.8038850177273362, 5.229142742704706, 5.150227939679976, 4.387917713549099, 4.089024405687361, 3.3168810813913483, 2.622525476712207, 3.2759705798060246, 2.153246346222959, 2.011097996145825, 2.1411743184367045, 2.1198690170693064, 1.3396966654015707, 1.173662537473211, 1.2368012851932209, 1.4918421916904463, 0.8180398138861136, 0.5410468264413889, 0.0]),
array([9.122193529260405, 8.24241566520224, 8.957761402884723, 9.088637478701473, 8.496861796863165, 8.48934897405461, 8.059077840954242, 8.017677163484942, 8.620142773805865, 8.572445566790595, 7.241996002217717, 7.024142807695581, 7.185670531166023]),
array([8.28605392652749, 7.763038743221727, 7.487871500424828, 8.12319159906873, 8.148432202269786, 8.06323997338348, 7.654322066101544, 7.719466327642752, 8.329954125791001, 8.085739383996522, 5.896883978591635, 6.909551997886716, 7.002195922419101, 6.846230718623694, 5.396744862983873, 5.614639154282148, 5.665074814790313, 5.4703265209868785, 6.796544377342117, 6.3356605621201405, 6.11123088384772, 6.491362231213037, 6.1408221845382664, 6.067273859814076, 5.177721166718934, 4.767636153079963, 3.7038523268647303, 4.75413540599134, 4.672246866715647, 5.2665708625456915, 4.100948955035584, 3.735286932825046, 5.140154637511036, 3.7104819213860973, 3.7588594797235584, 4.107723356566153, 3.0975744590784027, 3.018634026623602, 2.6567933715491527, 2.8101369545598747, 3.49715654398381, 2.8611646375685416, 3.583462234546955, 3.009782481715169, 2.9130223826487263, 3.019795026959124, 1.9888170721899778, 1.4296492778928012, 1.588865767946071, 1.6512769160265934, 1.623055656804549]),
array([7.609224356253485, 6.439732710590841, 6.884505613308216]),
array([7.316994178299395, 6.45656573449717, 5.952232845129023, 6.301324968843878, 6.244766569692154, 5.3732215436451725, 5.610075398166085, 3.3214440811686963, 2.7789879434845055]),
array([7.746391918220212, 7.768721024416119, 8.03752778016685, 5.631777162567803, 5.723082023025091, 6.464467793329213, 6.496700371450622, 7.218016400282018, 5.46511793600515, 5.498188768960562, 7.070729732971271, 5.683182798676838, 6.388099778281898, 4.4923934049442025, 4.852270791044179, 5.158778605469491, 4.0763969611960125, 2.9917311602510384, 3.072970179328543, 2.5145230504912064, 0.8193607555108606, 1.6560798265876744, 1.0124606992961311, 0.0]),
array([6.544208102481476, 6.011939109865599, 5.360296048936201, 6.458345513461847, 6.905195084795833, 7.240660050812669, 4.131329806472372, 3.9652045440406414, 4.9170938396617885, 5.044228775867661]),
array([7.594179638992498, 6.178953417286123, 6.930603007302507, 6.824642853085539, 7.085056581585491, 7.168710540611658, 6.344391875078593, 7.014918795252236]),
array([7.34887034459049, 6.671156382983143, 5.69308734366465, 5.575417025420055, 5.57179084900727, 5.629691922499993, 7.112419383989292, 7.210497133315291, 6.82054811339195, 4.515183894817876, 4.756570033385091, 2.9959800828097247, 2.598739698801542, 2.842197899264372, 2.967910239435154]),
array([6.858428055790584, 7.230590043338468, 5.881748611308718, 6.726050759897112, 5.638534378123877, 5.923901998946069, 5.330075647222509, 7.187534072088077, 5.43427323217249, 5.476806215660352, 5.693116499337374, 6.827750065267004, 5.110237881883698, 4.814222384552047, 3.8360542170166942, 4.6595385698720895, 4.758747236684735, 4.785701230149439, 5.263017943804828, 3.753258638091574, 4.908838583203645, 4.573459281216061, 4.325517710005889, 3.642678511191514, 3.8051004123988896, 3.981455959215895, 3.1263821193278476, 3.457857190078112, 2.666709382704003, 3.4032053329571688, 3.468662105168111, 3.5723405682622156, 3.4009158469048373, 1.8269959995852947, 1.3857708889668647, 0.7820191458867192]),
array([6.449074374315728, 5.506550335267505, 6.021659121234002, 6.65582039717887, 5.511523441990767, 7.018738536056984, 6.093776025168196, 5.6313693020656155, 5.560692587064812, 6.883107871766423, 4.813771787311398, 5.1488311808578375, 3.618334518554437, 5.1226646478603755, 4.0442564118232145, 4.430940301380438, 4.008484880078556, 4.771068861268653, 3.7784887607283744, 3.616753837154752, 3.998416212650443, 3.6653853506339886, 4.415465376096716, 3.15723025779707, 2.9256029731977296, 3.5937847930607254, 1.8824626088847038]),
array([5.838465850820517, 6.58033517702274]),
array([6.532834918254114, 6.856520007108268, 6.973479982122609]),
array([5.852426604722626, 4.5587730286324195, 4.663435188383708, 4.929609301340311, 5.189424710765424, 3.214849978462691, 3.053473036973331, 2.638132134970769, 1.5094062409606266, 0.0]),
array([6.082788087146426, 6.466047099222721]),
array([6.444354755709291, 5.463554268994628, 5.817420480478303, 6.187471305241228, 3.7420306796311875, 3.6715782533835926, 3.662147340348657, 2.828758983855523, 2.8491535312486302, 1.745434181608875, 0.0]),
array([5.695814466445663, 5.443386475687526, 5.58697356036914, 5.59579293276633, 4.475396305163702, 4.691704460877783, 4.665880916878006, 4.753958388364589, 4.176404532457892, 4.272876064207994]),
array([5.127088014739415, 5.026310400785705, 3.1719015651750455, 1.4757996302467449, 1.6731817113331402, 0.0]),
array([5.417995702129869, 5.5803369358565, 5.508365644740881, 4.495742824401907, 5.002954102425078, 4.16157602968164, 4.042159099833483, 5.255762873837186, 4.655732117063184, 3.8626260711711655, 5.150005995133362, 3.5710431785878547, 3.415956593756929, 3.0071458465698733, 3.2482415493023624, 3.152221275922136, 2.896077567106567, 1.9070050216899772, 2.398731309713654, 2.3644522335821883, 1.636753229586775, 1.5870355845244561, 1.5463793583330956, 1.020901348535964, 0.9025332430573859]),
array([5.388045747604206, 5.448971537236321, 4.499630121686186, 4.843554440650693, 4.700268346927876, 5.0024986391496835, 4.875711820107803, 4.433728772542657, 4.6924750429970175, 5.028436589109984, 4.945700084454287, 4.6039817256000966]),
array([5.173726603962839]),
array([4.448992740350803, 4.427842057396005, 4.0738207506485535, 4.597904070484591]),
array([5.0797804467099335, 3.3725851421284796, 1.520994061674611, 1.1099791877926708, 1.6459753748645765]),
array([4.754125536878346, 2.9187839136472022, 2.768676798923619, 3.4656741803125017, 3.5704824825796586, 3.575628779836782, 2.840836182089322, 3.140422671936033, 2.2632578305265714, 1.4090747438888764, 1.045749007739597, 1.109720211427622, 1.5001488997796746]),
array([4.502642705280227, 5.073690146431603, 4.536097726060506]),
array([4.872952849367805, 3.797790263051258, 4.741380865549605, 3.715455492814077, 4.51661965550985, 3.0680126416902467, 3.221737485486075, 2.6767428355597573, 3.1068002284360636, 3.450977261635683, 2.4726278768820147, 1.2367309226040781, 1.1291502781817018, 1.0188716209203426, 0.7908394429720234, 1.096757467362249, 1.247324229047651]),
array([3.4514011367880877, 1.5733093233445103, 0.0]),
array([3.7004342092093427, 3.846555447223752, 4.1625651222510225, 4.724258665534068, 4.2795071521662855, 2.675957968237978, 3.4314744909599115, 3.23935459934143, 3.330517012931229, 1.9170376329995658, 0.8076253499179697, 1.269953926293103, 1.0244468687925214, 0.0]),
array([3.955353014806236, 4.2109428666496775, 3.6532912239235262, 4.474335848931338, 3.7065641691827036, 4.5670208051679335, 3.946341224124402, 4.752456525755693, 4.3910204062197575, 2.738123694574747, 3.2324367356786134, 3.1667351964367807, 2.9926986945806, 2.218648940204468, 2.0877830004937326, 0.8348945269539663, 1.052123494485802, 1.7600011962438493, 1.0074205973167678, 1.7618679388312923, 1.5564808773098795, 1.7515482994769498, 0.3336778325777125, 0.07362262097031548, 0.0]),
array([4.678459363113648, 4.588024311272456, 4.527490181418553, 4.493408029058543]),
array([4.468843324236254, 4.518708789678376]),
array([4.172759585654287, 3.934246548065205, 3.905281547434685, 2.91486938939734, 3.0360195263527467, 0.0]),
array([4.415269789073326, 4.347455818908384, 4.463807130040287, 3.996107897362526, 3.907082325043974, 3.8885186860091814, 4.3681345722796605, 4.436446322045786, 4.01819016482391, 3.2645677692820163, 2.9579505136426767, 2.6707634092665447, 2.1919868849205884, 2.4069576692470553, 1.1220825745460759, 0.990604628311154, 1.0604097248452007, 1.1199518155806465, 1.6354393144058792]),
array([2.711897472240282, 2.807145628855955, 0.12232167469055379, 0.0]),
array([4.079327895797182, 3.664666959601604, 4.295131839406057, 3.3951876119821516]),
array([4.304241751703268, 3.8165626887166963]),
array([4.169880194722375]),
array([3.912051201955532, 3.9798981106057134, 4.011267456123417, 2.3911817839985403, 0.2291126392122087, 0.0]),
array([3.5246151249276885, 3.5794662070491747, 3.3880471238994048, 3.0091216526538336, 3.429136602840896, 3.2923906999293564]),
array([3.9844419660512047]),
array([3.6880378181358826, 3.59283161513031, 3.4459305407781455]),
array([3.7676488456553585, 0.03170400852842167, 0.0]),
array([2.868986669401189, 2.8709552785809205, 2.640657847462328, 2.554058565701454, 1.028548664208869, 1.7869780854722375, 1.2445772626599714, 1.0199838475179983, 1.0181467089269687]),
array([3.272325760859842]),
array([3.0889793389382034]),
array([3.5053096392726837, 2.2976299671857783, 1.4317354468398784, 0.9560495032816515, 1.7865652467212416, 0.0]),
array([3.463594569918932, 3.2761593297528857, 2.7429743887730935, 2.672868663991993]),
array([3.3197287766191432]),
array([2.9301628759366034, 2.357032907356976, 1.2956315065335229]),
array([3.052128549743184, 2.852083948717557]),
array([2.877418217900665]),
array([1.7917681846870661, 0.3970122220230113, 0.12395566421976248, 0.0]),
array([1.627931870436281, 0.0]),
array([1.584728533098585, 0.9755781188393247, 1.473499717814827, 1.4755608599228072, 0.7985598444383775, 1.1384409745639452, 0.892556552528591, 1.5445150413563837, 0.6152485738639046, 0.0]),
array([2.3674118137568536, 2.392466335285082, 1.9172032697552817, 0.880359505185512, 1.0245826265578568, 1.5286626711392155, 1.0186024547830024, 1.6530812117056428, 1.5569865880547114, 0.9904935924415074, 0.688203837058437, 0.0]),
array([2.0890749965313296, 1.7805591297506358, 1.517180125052687, 0.8964215752263086, 0.0]),
array([1.4315205116959673, 1.1336375997233006, 1.4027904027438989, 0.0]),
array([0.9508605847627547, 0.6393683980946376, 0.0]),
array([1.2695563811050063, 1.2608408766086026, 1.1003079304337757, 0.8150116859376864, 0.9517138224178444, 1.3646042319699974, 0.0]),
array([2.0411576080706646, 1.4164836133326022, 1.6346173106758608, 1.3250052912992487, 1.5628481752215782, 1.5988135167493664, 0.7679767855354321, 0.0]),
array([1.7600685136567542, 1.574813577186627, 1.7571215408491339, 0.8058051456867292]),
array([1.6159740235023103]),
array([1.8372884426260716, 1.5846754577409738, 1.22351539261029, 0.9605975527955529, 1.6357122079796262, 0.8870787783227075, 0.9811303937606178, 1.207013513716392, 0.21389208209523036, 0.0]),
array([1.6769335838056527, 1.7585027772816997, 1.5954995237976968, 1.2803432581728353, 1.0208383885425416, 1.4635490521216727, 1.2336751364365615]),
array([1.4466575982457004, 1.665528873393539, 1.1475371281501994]),
array([1.1289425119464835, 1.5711815462116796, 1.0589528241798294, 0.0]),
array([1.405284882048368, 1.477200417856566, 1.4977866540944884]),
array([1.5405318579640979, 1.6656518997873218, 1.2706457657788148, 0.7710851240795927]),
array([1.5410858038051638, 0.8820019460566912, 0.798244148758321, 0.4976447326386477, 0.05328655144011356, 0.02554929742261798, 0.0]),
array([0.9537443662755489, 0.0]),
array([1.5506075888155362]),
array([0.5013539673218056, 0.0]),
array([0.8627372320462597, 1.1834088917311312, 1.148784876000508, 0.7389580738770813, 0.0]),
array([1.2130508277879881, 0.0]),
array([1.172891334029932, 0.7706987049410862, 0.0]),
array([0.4317613460114923, 0.0]),
array([1.2666221935344903, 1.428217323771072, 1.1629521911767973, 0.9075280039752711, 0.9280501844876622, 0.8241231126148459, 1.0716885439440231, 0.0]),
array([1.0517276594884724, 0.0]),
array([1.0653214625462555, 0.8120977411170607, 0.0]),
array([0.877878449896806, 0.0]),
array([0.8070090032872268, 1.0883335710822244, 0.9604798084619703, 0.0]),
array([0.2090701384738063]),
array([0.9727880353880007, 0.12808429976961366, 0.0]),
array([0.018233984617510093, 0.0]),
array([0.7220344457943326, 0.09459474654752525, 0.09632315504199171, 0.05761804585860629, 0.0]),
array([0.8382830760212476, 0.8932911014609621, 0.9152286277221218, 0.8779335908321756, 0.7378141518950366]),
array([0.7900821958376503, 0.0]),
array([0.4536377457782243, 0.0]),
array([0.0696939596868564, 0.0]),
array([0.05618025212188532, 0.0]),
array([0.15968666419170624, 0.1752358216678363, 0.0]),
array([0.07100223458204136, 0.0]),
array([0.022920797340432925, 0.061756867635585716, 0.0])
]
d = [data_1]
names = ["29"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T13', 'T14', 'T15', 'T16', 'T18', 'T20', 'T21', 'T22', 'T24', 'T25', 'T27', 'T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T35', 'T36', 'T37', 'T40', 'T41', 'T43', 'T44', 'T45', 'T47', 'T48', 'T49', 'T50', 'T51', 'T52', 'T53', 'T54', 'T55', 'T56', 'T57', 'T58', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T71', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T83', 'T84', 'T85', 'T86', 'T87', 'T88', 'T89', 'T90', 'T91', 'T92', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'T111', 'T112', 'T113', 'T114', 'T115', 'T116', 'T117', 'T118', 'T119', 'T120', 'T121', 'T122', 'T123', 'T125', 'T128', 'T129', 'T130', 'T131', 'T132', 'T133', 'T134', 'T135', 'T137', 'T138', 'T140', 'T141', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T149', 'T150', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'T158', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T170', 'T171', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'T180', 'T181', 'T182', 'T184', 'T185', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'T193', 'T194', 'T195', 'T196', 'T197', 'T198', 'T199', 'T200', 'T201', 'T202', 'T203', 'T204', 'T205', 'T206', 'T207', 'T209', 'T210', 'T211', 'T213', 'T214', 'T215', 'T216', 'T218', 'T219', 'T220', 'T223', 'T224', 'T232', 'T233']
def get_taxa_names(): return taxa_names