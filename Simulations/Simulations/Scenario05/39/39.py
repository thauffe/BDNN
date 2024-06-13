#!/usr/bin/env python
from numpy import *
data_1 = [
array([31.181397962435227, 32.439377311310956, 32.6487898737647, 31.954098341680055, 31.074134522210667, 31.880530988420606, 28.568158921178238, 29.48407273793374, 29.838107595069026, 28.226668710743745, 29.55843532004132, 33.65049944533641, 29.012936068255357, 33.696650872986865, 30.179404358547643, 30.835056766287792, 28.2183507259782, 27.145656874923397, 23.265878447776892, 28.051220136194658, 24.28718414174888, 24.246589190951205, 25.80892044905613, 25.956726674876265, 27.40916629580078, 23.695987830453486, 27.314987240530755, 24.22899728055506, 20.865034737655673, 21.851250212055433, 20.59608664712743, 20.98529619272345, 22.095394345044483, 22.900785446978045, 20.7490076687998, 21.320909734000917, 21.190512965135266, 21.79860255775019, 21.612399203620114, 21.559978476534823, 20.696899910963243, 22.53104634719626, 19.0597538995869, 18.084515289540146, 17.4052712828724, 16.893815203565737, 19.11483957361479, 14.647441138260634, 15.717137546459615, 15.24697280247892, 14.644570980315063]),
array([24.558849459646037, 25.604939724475713, 25.211083056559737, 23.36964055983471, 24.585699853364552, 26.27682943536703, 24.298311781331815, 25.77709623896812, 26.583134068872866, 26.33404689765803, 23.89409031370839, 25.376797725717676, 24.530070383807633, 24.28512667553946]),
array([24.690731813187217, 24.92511191351783, 24.8148170176538, 23.64354173555914, 23.870844169234864, 23.589818994139073, 23.948334993941362, 23.746773142355607, 25.09895400672102, 25.130055304886785, 24.461050436700777, 25.513027185588257, 22.72896436060923, 22.43624722313842, 22.40086042408241]),
array([23.460366751247914, 23.142942713933607, 24.017846664156608, 24.03707694812241, 23.650381275219058, 24.678569401486243, 23.859565347579885, 23.514020623737068, 21.6103491237299, 21.442940884503937, 20.574048187434748, 21.42882790748122, 22.639998925886644, 20.99619574916413, 21.077223771583377, 22.007590423939806, 21.790342974676463, 21.96153572202491, 21.345549696419273, 21.058879926733773, 22.09488107025257, 20.87206176485011, 22.821784044565938, 16.93379415706694, 18.80766457770212, 19.707502605980636, 15.874947341294616, 14.575825696118695, 14.187489371135404, 15.826438250313059, 14.511964311141746, 11.825260610100628, 13.147813418082189, 13.526271879392084, 12.633030965085423, 12.279521773176318, 12.983705946506864, 11.962329071476773, 10.280506674181632, 10.165983857185841, 9.81890455851241, 8.239851002722741, 9.89770006288064, 9.682834514895212]),
array([24.17172646808297, 23.66515193223408, 23.762646013021342, 23.427718431107646, 22.62612969800464, 21.719865330099964, 22.23634930743383, 22.32487320069258, 21.582903404239282, 22.80024928849215, 21.27455363275163]),
array([24.145271847576442, 23.823967609854023, 23.859402745065157, 23.47636045314217, 24.18159801983172, 23.108751263910747, 23.615611281750954, 23.750302089922847, 21.73883049950447, 20.74316783389284, 22.72951825542802, 22.529677382120113, 22.943553367335234, 21.37719749396382, 20.99633378503393, 22.38835731149528, 22.830440316783104, 21.994110080668268, 21.03353588614965, 20.79986776321019, 22.119763114116136, 18.621439133471746, 19.016149789886416, 17.484836278203545, 16.701865633716068, 19.908414467197147]),
array([23.491904529676226, 23.157964714414, 23.041380553200913, 21.89418151028757, 21.49698232198026, 22.56700317344914, 21.881775275506964, 21.950725282448687, 21.61883224416857, 22.42734406584972]),
array([20.905227351138358, 21.758925938858123, 19.607288725650875, 16.903330928391753]),
array([21.90542142928417, 21.033203168590653, 20.93687708856688, 22.376303178212947, 22.122106628472153, 21.126806072321315, 22.514618976318445, 22.32930933277983, 20.45773575475796, 22.398627602987027, 21.44552574406769, 22.2704944295834, 21.04730599163719, 22.28095663353904]),
array([21.89449746811867, 20.779804020828376, 21.69879867179519, 18.448959387964628, 17.590581300939352, 19.82242435189826]),
array([21.411471974514473, 21.430528006275875, 21.644807953298447, 21.649893778793402, 21.84198871490563, 21.79612815831973, 21.662791041775595, 21.895005674120284, 21.69098583387849, 18.417495902306, 17.620735825166275, 17.341874944041678]),
array([20.981619019085883, 18.809988426328957]),
array([21.251065686514963, 21.45276993599324, 21.077470066049802, 21.46841161393145, 21.515067691134952, 21.028250544432726, 20.856397023775873]),
array([21.404873750973678, 20.69697952431166, 18.089799424851467, 18.06288040292195, 15.727167273302433, 14.907302221323722, 15.71895773922396, 12.827343378640668, 12.468390511922514, 9.890661192498458, 7.991936668665909, 9.016787719585807, 9.902552020179542, 4.856034681214467, 4.706117432973191, 3.055848285675953, 3.0063311066528366, 3.065944880081912, 2.8281478388910974, 2.4010760480231292, 1.1129104509318317, 0.020966382327541444, 0.0]),
array([20.652191469908807, 20.971485745467586, 21.040305173500318, 20.184721548467294]),
array([20.07167433088863, 19.875136122303534, 19.1173595319925, 19.08324849007971, 19.978649426263207]),
array([18.474949039345447, 18.285962643528627, 17.800508938790635]),
array([13.89607485473719, 14.751092701147517, 14.040674451021303, 13.300639617901853, 12.913623476047624]),
array([18.56997592138102, 15.070812925201555, 14.966278246719506, 15.567759228400748]),
array([17.561839021619672, 15.998177139394167, 18.23703518003774, 15.065828786611306, 14.72225258014133, 14.72487564470883]),
array([17.653850475003303, 15.821978386957948, 13.214665737909076, 12.081981996980764, 13.48988244270934, 9.53232686229969, 10.79621170741809, 7.730073524891603, 10.892963345246018, 8.960470492863088, 8.588539658120466, 10.934741591955794]),
array([14.224052403410012, 14.946402111409975, 14.276090900713076]),
array([17.778468971286603, 15.833989324290558, 15.643104620602442, 15.747014511447432]),
array([16.26563092801648, 14.80414667854617, 15.789827998205245, 14.141861783840396, 15.278223000764406, 15.026745644457426, 14.313496664690703, 15.633641230071923, 15.821070653375212, 15.550702492795114, 15.324336552025748, 13.072064434600168, 12.956682359651657, 11.805035971358368, 13.680657070228534, 11.75590546469249, 13.400713803841231, 12.659635442189211, 10.8572145928099, 9.871016989915901, 10.797140387289932, 11.506685140383706, 8.287626003826539, 10.416359434495627, 7.445553593195959, 7.564400317587941, 10.462053447856563, 10.540159290818732, 9.364109694017191, 8.740695612879355, 5.337734021744334, 6.931838018728104, 6.337289014553147, 6.680354490363021, 5.678472094894657, 5.943220135784961, 4.4450385569008946, 3.7096628418811184, 3.7861767095529055, 4.651265346126046, 3.3866741910608322, 3.2160447888982895, 3.307278439502402, 2.662862994000377, 3.225705129697637, 3.1464451875581387, 2.0146556644529814, 2.4683058007011613, 2.176766729976826, 2.1408386134990964, 2.446485015449437, 0.9992110291648498, 1.281990459597195, 0.38865379072822187, 0.16663806479158005, 0.5794708664494693, 0.0610727979868943, 0.0]),
array([17.604821725776738, 16.90029244812888]),
array([16.185483485139944]),
array([16.90422007948711, 14.861762796416897, 12.465059062197849, 12.433785865069769, 13.108811484195217, 9.450688295102555, 8.971160428722474, 8.43658811669939, 8.016235662264961, 6.169743439816381, 6.370770932257012, 4.956752964465252, 4.063805134588252, 4.735781703590305, 4.614577962550951, 3.0242950800455146, 3.274983047257247, 3.2991464487002102, 2.12921043180402, 2.5288799069837133, 2.182135021615938, 0.9325617264723356, 1.0276989085966024, 1.22813801320842, 1.1172841259060369, 1.0087998927028328, 1.4075120984580534, 0.7824095344102584, 0.0]),
array([16.918480785096985, 15.499161579335299]),
array([16.774338965622196, 17.060645841087354, 15.530559904976364, 15.37819680161458, 14.568946075905925, 15.627636821436344, 13.80742543060018, 13.569140878707088, 13.529148322154027, 12.898827039508543, 13.499252020676645]),
array([16.55027750682632, 15.371480692000064, 15.489078378637869, 15.767060841962381, 15.709851892664508, 14.487052533928672]),
array([13.337906309707071, 13.434092611534172]),
array([16.32825433241385, 14.875440089418456, 14.487384556880995]),
array([16.3210340911241, 15.68399166285498, 14.82246417901912, 14.356734771727718, 13.733214674812533]),
array([16.03215910398176, 15.014523309389249, 13.862906990190478, 14.45388245960973, 14.8454569293226, 14.372054639135836, 14.484694954956682, 15.055181945314942, 15.752991333471929, 13.91907965994589, 15.83369005432112, 15.727744092062922, 13.244955099639911, 12.763791448474723, 10.287869828630075, 10.82152843167118, 11.071014892716228, 10.174615135555378, 10.477669615972879, 9.909040206651047]),
array([15.9605309832726, 15.353587863386457, 15.101697982347742, 12.81695506289592, 13.729620699708676, 11.698360398255497, 10.374634264158876]),
array([14.936913635195149, 15.420869101241179, 14.394551313441127, 14.109104823475946, 14.826484565233184, 14.04604580989093, 14.574171479293172, 15.743313188040567, 14.345268041230096, 13.889185385156809, 14.290270058869588, 15.936642997816133, 13.592161590616493]),
array([14.856290860305203, 15.042798584618852, 14.198732976878226, 15.469033510710878]),
array([15.233606251415022, 15.65710286525584, 15.687078495988636, 14.632947455189461, 15.650560098050981]),
array([14.781769678280037, 14.43429177220991, 14.104126797915468, 15.639912493678056, 12.375626401010006, 13.139321211242876, 11.361535276810294, 11.507500037732774]),
array([14.637568733479704, 12.038648029904454, 13.28305415463554, 13.416989203598892, 13.532155485705795]),
array([15.022388310456567, 14.98460792732167, 14.237009744527702, 13.929660959527231, 13.18249385006266, 11.98483268814111, 13.195037068545542, 12.839149192599974, 13.684881695695376, 10.054381023168695, 10.113417128896657, 10.481771343780192]),
array([14.860083916425634, 14.93830290710688, 14.943005093143455, 15.026491166962984, 13.550891197864388, 13.268681407547248, 12.881907652416963, 11.678213676036197]),
array([12.286107344362948, 9.028509422054455, 10.922864797644811, 5.980058772334547, 5.74967485786885]),
array([14.773878461909943]),
array([14.717234068033807]),
array([14.298544427911398, 14.38013316136863, 14.147092562464485]),
array([14.411105766394911]),
array([12.919147962366916, 11.83667885331463, 10.766805567480286, 10.95433323192691, 6.593495227310318, 6.685166105879906]),
array([11.168512020734168, 3.192356006765318, 2.4311443426088437, 0.9848241257925328, 0.0]),
array([13.551982820858449]),
array([14.15206713098764]),
array([14.451509492081884, 14.24269865930104, 14.238904765907945]),
array([12.3677911000892, 12.195612847903568, 12.936452763692095]),
array([14.155279082969011, 13.965446222261082, 13.705036683909313, 13.010559476259646, 12.467819141691988, 12.621542835303158, 13.482548184768879, 12.772646596229434, 12.136332436452284, 13.41476472139803, 12.507996047636494, 11.588612706672608, 11.563392747060224, 11.54988035061458, 11.384171474822582, 11.57742895739342]),
array([13.912164871685759, 14.005708733928728, 12.1733694206298, 13.526132726668509, 12.408195553617501, 13.442061944359137, 11.970885610513282, 13.328714881074447, 13.250878663119165]),
array([12.801171546824207, 13.131373460386257, 10.805154323750191, 10.38950360398037, 6.526313191729809, 4.414205193160837, 2.7207024030905975, 3.058342251766513, 2.042233945563098, 2.0083573722726102, 1.3001623404422735, 0.7108644353211362, 0.6364767438881659, 0.3755411375270069, 0.0]),
array([12.497006878615949, 12.576284727851249, 8.506562328232363, 10.891678747912794, 9.298257999841761, 10.59422778609368, 9.835348095051359, 10.53211767372408, 7.696834278193769, 4.3690068097921975, 3.5154117223330794, 3.568596516785285]),
array([12.578976990997726, 12.82347108971966, 12.557515867094875, 12.782769170728113, 12.11755207218512, 13.620975213687114, 11.939222267593104, 11.300736981520124, 11.507326637666212, 10.415883073631843]),
array([12.570111590969933, 12.711360044899523, 12.993828231098227, 10.458632508783447]),
array([12.95168983373395, 12.652150778870878, 12.14180915198108, 13.325456135545403, 10.957203362490926, 10.783439522322768, 9.603534605443492, 9.850478057477174, 11.628410259148529, 10.948893429145459, 8.240626223002568]),
array([13.225951838030424, 13.05679489363321]),
array([12.699482483032963, 13.255769368386892, 13.20714682947042, 12.873427269059526, 11.849015997004468, 13.476273991803186, 13.522335194592813, 9.632625163341048, 8.062525089631384, 7.710102715346477, 8.304541251456486, 9.180655579669093, 9.759469464392069, 7.502964904182262, 11.463419504980358, 10.663872887180377, 8.6260810492229, 10.750236950927427, 7.828062071803293, 10.107148892130937, 11.122493777963193, 8.863754081321018, 10.930583767466135, 11.442205351265589, 10.040140344326373, 5.621602080330349, 6.511968698077339, 6.9051008007381665, 5.359660510061953, 5.358904474372385, 4.442904737858319, 3.9285586438046787, 3.819415054637058, 4.996921657793872, 3.892130360350876, 4.033362836691509, 4.8829551240167435, 4.560610620723946, 5.107492506481525, 3.1045463377703157, 3.199730288797058, 3.411518724361968, 3.3862604992914522, 3.085640215400085, 2.985856541562223, 3.2960431088065048, 2.789291581790372, 3.231748158258207, 3.2569940418792975, 1.8475384383814868, 2.570782626395853, 2.4941645699207893, 1.1619282568899014, 1.4186832674850267, 1.012513557862937, 1.058612491657763, 1.3323646269012854, 1.1170259653396422, 1.2096052773136416, 0.8154866041431851, 1.1918191393632565, 1.6335307231050962, 1.619163526164185, 0.44942377026484354, 0.6596285916584903, 0.13639826129024202, 0.21340666791201135, 0.24686821982180107, 0.0]),
array([11.864543056341992, 12.877805572921904]),
array([12.510742553878304, 12.800723748936202, 12.021345965102322, 10.82167769394913, 10.443453194646818, 10.426389618735332, 10.543004158942141, 10.64816647582646, 7.619901637855717, 11.10526368831027, 7.5659812815071605, 6.136654663996283, 7.074455610892529, 6.269184084151771, 3.940677979118017, 5.004552810992417, 3.975464688162567, 4.678153971940864, 4.470824028458381, 4.650173559215036, 3.1499665497106046, 2.6068142650071153, 3.125124395052409, 2.6934771562361854, 2.84484898269589, 1.8909478522168757, 2.092717323736916, 1.97140998893515, 2.108410848460528, 2.1691550046465715, 2.233798795018603, 1.5868218352881374, 1.3129264511401588, 1.070638360803818, 1.180197783162365, 1.7403042811419813, 1.0931004164759082, 0.21390020103615714, 0.0]),
array([13.417766386558773, 12.93093348466977, 12.944825935819335, 9.673617994440027, 10.66120379484131, 11.24899845290331, 11.231121738206452, 9.200819564658257, 9.523192011228215, 7.026949573883225]),
array([13.232286800161663]),
array([12.664469945719297, 12.259756165868994, 10.834028506714732, 10.446763145506754, 10.420362604147588, 10.89753858051302, 10.360308025508179, 9.693010797196047, 11.520312126487072, 9.28768896542232, 9.027474786326515]),
array([12.983692190182548, 12.748815686754902, 12.518114901914608, 12.822861431444759, 12.686842892876285]),
array([12.709090163841866, 12.861214134375553, 9.892038049923391]),
array([10.74910172227365, 9.833973890414029, 10.429204049194336, 10.781691720949059, 11.158280464759947, 10.40502164724412, 10.871341486955195]),
array([11.99763664932704, 12.160421235509835, 12.38694327064298, 12.745363963965097, 11.577321623797324, 11.594540021211191, 8.839968993464883, 11.194354218462022, 11.573028454581447, 9.531607444808387, 10.07153735684133, 9.893801492193631, 7.022310208648463, 5.833422408939369, 6.876366218741506, 4.79058158622429, 4.548344740154015, 4.118432121418517, 4.646338612960976, 3.172718731254466, 2.4386051643428135, 2.185938938674974, 1.719324462095829]),
array([12.861375143641258]),
array([11.998060480054702, 12.637340951261447, 9.704600546519567, 9.188665762246963]),
array([12.067004194844898, 11.974285310907936, 12.389005529040643, 11.740257696487292, 12.384308827011166, 12.53361102884061, 10.484199249177347, 11.340154139516516, 8.174630385019501, 11.333470940463362, 8.073133944930378, 7.865628344849595, 8.377990306009913, 9.538734705067405, 8.09769209444539, 9.59374472255684, 10.008292462304361, 8.747970931931144, 9.690194716642074, 8.864272735119187, 10.895355931574633, 11.415150717661989, 11.048872257958804, 9.436112568551476, 8.07952819012969, 7.600497553387156, 8.715636919821726, 6.296698501132056, 7.119097768558501]),
array([11.977302726470525, 11.469294638447979, 10.843700699250949]),
array([11.382783268322504]),
array([7.868210416461574, 9.205642081348486, 9.161420203911034, 7.288817095381348, 10.577720468282845, 10.741297098417217, 9.847820631732324, 7.582567821313376]),
array([9.774815785422188, 10.166093369447841, 6.609133646478742, 6.808270072416264, 4.047777771016584, 4.38174672020906, 4.273409877528791, 4.616693623747758, 3.4347042000932917, 2.7177396587723517, 2.730758610165408, 3.0814970854232984, 2.949136884707336, 2.0758276525923156, 1.4234875382882304, 0.8391021784801366]),
array([9.91791123910372, 10.85123750255246, 11.258335410055633, 10.855509515278706]),
array([11.930411240800472, 11.846638831285539, 10.530301809471709, 9.921119204181792, 10.660404808613873, 7.940023540934741, 9.941050333006594, 8.12597704704979, 7.5984451258412955, 10.535260527782716, 9.36810771109003, 6.963271992743084]),
array([11.929019478569113]),
array([11.770660033610925, 10.329608118918681, 10.15796060038229, 10.708255192159278, 11.393887528251863, 11.156680648969273, 11.411897533025039, 10.662442033759728]),
array([11.661658984831556, 11.149497476054513]),
array([10.263052974868085, 11.356476184918172, 9.570426373103025, 10.111584532158744, 8.541523006978462, 11.330498058105832]),
array([10.94632367703998]),
array([8.220796932817409, 8.7945787364421, 8.512507077252323, 9.088899024695356, 11.079933739178596, 7.989621933422192, 9.978745338761803, 9.420638555894985, 10.522542789487465, 10.5391324875099, 8.491906264412243, 9.91983334633981, 5.74087994584624, 6.965153891369427, 6.16819178113393, 5.831692122241751, 6.6329125180638355, 3.69056830638071, 3.9839497009236267, 3.9035471380628857, 3.9890538154683735, 3.3224842821549045, 2.7806196224533535, 2.7903938808409716, 2.6064640379403814, 2.91201804965788, 3.5188576081231364, 2.6492234560001235, 3.584129178077762, 2.0875025918168393, 1.9185691133202007, 2.404455334941548, 1.8220524340145081, 2.10596184613495, 2.207332667319624, 1.5670885512699473, 1.0676668869043713, 0.809915412713288, 1.469717194305461, 1.0064175311788892, 0.8528218712120827, 1.6912847388395282, 1.7784769942964216, 0.8757799351273254, 1.5792516667531762, 0.46812100703951387, 0.6899241628998716, 0.5981120776091271, 0.0]),
array([11.384122636791325, 7.778184269716167, 10.534036916632703, 11.126164326091363]),
array([10.910985236698817]),
array([10.059123416823137, 10.354970732561684]),
array([10.944984944964842, 9.492217597598913, 11.297434120527381, 10.595530747094246]),
array([10.198119870442385, 10.589912220430636, 6.991991784661618, 5.33808884926176, 7.071261529924022, 6.343668938353484, 5.90845840513296, 7.018746118680479, 6.609649057395646, 4.00899216490707, 5.068661287811892, 3.735204337388951, 4.422126287439248, 4.294344657073045, 3.68944305471743, 3.8374627862601898]),
array([10.926414713604277, 11.113537205571918]),
array([10.05875306091809, 8.851376161193938, 10.088420584095173, 10.648095990212255, 10.59760954873486, 10.992784045615215]),
array([8.812683423111247, 7.816678773406606, 7.8442791845329065, 9.93517425271315]),
array([8.877339319467852, 9.651199959366274, 8.101526452109301, 8.056368439827276, 7.926138307489759, 9.069694215082285, 8.156481253737203, 7.033813990246453]),
array([9.405904280624744, 9.127366771558664, 8.912793768647731, 8.40599414200885, 5.687724101975046, 6.472330474016809, 5.185380534187058]),
array([9.006033126897629, 10.672909776216049, 10.213203300266725, 9.978358384025782]),
array([8.292228899045258, 7.845091707826274, 6.330190753652719, 6.198108608676864, 5.662811108677041, 7.161817887351711, 7.154783570382313, 3.390428080515539, 3.2980769785551827, 3.409538465090738, 1.9709446399547463, 0.8575910302763391, 1.3625759726285214, 0.9678073430575002, 1.1038435676221634, 1.0377489574551584, 0.8678227189664434, 0.0]),
array([10.127809977116938]),
array([8.646380761598525, 9.337194885517906, 8.559971185926926, 8.595744953146792, 9.674677702651426, 9.957968158114499, 8.612702640523986, 8.46434982119731, 10.276472232012418, 6.788354468925205, 6.2625362986434805, 6.811816303179987, 5.07556845534206, 4.968229345766572, 5.054459234811762, 4.847000312766025, 4.9328765312571035, 5.249286730729063]),
array([10.087598931012211, 10.32810641530298]),
array([9.997840373852796]),
array([3.0876530054761675, 2.1140279554192487, 1.4368702053817695, 0.0]),
array([7.588462540412552, 8.888737962572453, 7.441869272871519, 9.388046211926149, 9.913141480944859, 5.874974788537193]),
array([9.363517176370252]),
array([9.205031269535588, 8.62726530948205]),
array([8.995792884850704, 8.890356613588487]),
array([8.989854397839547, 8.874374396958357, 5.710986351306497, 6.215012042825118, 6.54501483620504, 3.9865223710310973, 4.95635692892475, 2.67477615437046, 3.2628791435779654, 2.660762688199425, 3.256015794592063, 2.4604970251073492, 1.9777932167539414, 2.327870209417017, 1.7189397054602256, 0.22113209362307307, 0.01466783376288762, 0.0]),
array([8.593419311468423, 8.570409132573916, 6.809429764861759, 6.088152647243806, 4.568153576542456, 3.5741743777119326, 0.8205211179130417, 1.2186425787935158, 1.556631154136254, 0.0]),
array([8.23604129890108, 8.041452398563141]),
array([8.353415654811036, 7.569810829882173, 8.72248281143084, 5.249703456086303, 4.385566476403221]),
array([7.389372470467527, 8.104526975124095, 7.2900913979639785, 9.311022452646457, 8.27889971788995, 9.012837405513787, 8.538600812683637, 5.938564642371867, 6.379726238759308, 6.994578384142334, 7.054727279545936, 5.626579164719694, 6.407735263116081, 3.8855468435546063, 4.110706497213147, 4.851295165422517, 3.733392944053402, 4.6574697972976065, 5.236491121163248, 4.556959127551584, 5.096717297952502, 4.8155952586210375, 3.791306972492053, 3.6937959673481076, 3.6809881364667136, 3.5663783725124705, 3.1851919270862896, 3.3054976840300228, 1.9582870378142343, 2.22487967803338, 2.1568322303675487, 2.321787026038005, 1.966577636138116, 1.9003771773478286, 0.8279152214277802, 0.9521184419562544, 1.627131117756654, 1.694824049092873, 0.8672735989381107, 1.5265872756263632, 1.616937627353598, 0.992896579655247, 0.7326053233249694, 0.6235969216651309, 0.0]),
array([8.581069594083646, 7.870011421470474, 7.004426221632076, 5.71307733294238, 5.618709681613053, 6.433136880009382, 6.891057296935536, 3.759177509574656, 4.658548372963072, 3.9883109242569503, 4.920614298556285, 4.55178528570655, 5.283378631180445, 3.9235889905324783, 3.8249253538265546, 3.271691698755448, 2.5849924578070063, 2.7537680282779204, 2.8960218454680193, 2.9197530482590786, 1.87572863334736, 2.194932821967765, 2.322069121148509, 1.3616063618334382, 0.7972025225790611, 0.8768263238891341, 1.3608731678354555, 1.1851927036041823, 0.5387140653033773, 0.12333491830307523, 0.04466116499825827, 0.0]),
array([8.31613861129636, 8.642855939687706, 8.050382196942419, 7.015521580734321]),
array([8.684135477953573, 6.098822547646216, 6.553379056793942, 6.012889751465622]),
array([5.829856989378928, 3.6782098589227017, 4.597729435630751, 5.24227166622213, 2.8422053720965055, 2.1880036128788714, 0.8895573879246674, 0.0]),
array([8.454941716317983, 8.427219050580772, 6.1644997559105255, 6.746391147509523, 6.45206939891461]),
array([8.3291944454166, 7.7288967665610055, 6.627863584570064, 6.384701760287467]),
array([8.278379716896875, 8.050676907746823]),
array([8.22819813149903, 6.987083591962275, 5.867013260025013, 6.020667350386507, 6.713601078156169, 6.35947425721228, 6.681235816526562]),
array([7.999115651957142, 7.7508421190236145, 7.809798308589912, 8.44710079811709, 6.773643071155001, 6.27290750147165, 5.639013310487729, 5.580382068137448, 6.995905679978399, 5.777106186418034, 5.561766947075061, 4.38825153651494, 3.9542548793998664, 4.76496918383555, 3.7145389739156998, 3.7458849609905194, 3.602243400849914, 4.401427642138334, 4.8936391901894005, 4.221195652600897, 5.2907184833300205, 4.045369603998896, 2.6734320599137353, 2.708088529390492, 2.783153839897997, 3.0003139284575076, 3.5661630283325665, 3.5351451345513714, 2.842771778171688, 3.4138065503622235, 3.2674856989541943, 2.934906324118855, 2.8421484264090733, 2.8265276489581295, 2.4215769298768524, 2.44502929708556, 2.2959955827187124, 2.191781165597979, 2.399457089059747, 2.2187574784113875, 2.316105291142902, 1.1010451086750273, 1.1958394523954716, 1.5575722405630557, 1.3584852003532069, 1.4590581509019043, 0.9204641469268967, 1.6444777490532823, 1.464981669245724, 1.0884532842863757, 1.2853093953382795, 0.19805733881753362, 0.5778287878006804, 0.37597760861988855, 0.21684112328483107, 0.7426984079525574, 0.11763482395606088, 0.032388573781910734, 0.0]),
array([7.308919701399594, 7.726598941997562, 7.118980553233858, 6.932460438949279, 6.810350248766587]),
array([8.147261160286801, 7.926273372006755, 5.6045441522588675, 6.474319780605357, 6.686923837672307, 5.971425772231719, 6.042426230498879, 6.054806264146327, 6.8451997442826595, 5.301382539839298, 4.675732493512459, 3.953196565638719, 4.977375615106801, 5.233630132041688, 4.455883192519995, 4.467299161458547, 3.954447385206956]),
array([4.6450203909083605, 4.470215652292854, 4.886540263950383, 4.138652023558128, 4.334636346520711, 4.864167105914306, 4.19746061405659, 3.031846097469229, 2.710146624302248, 2.806136851628013]),
array([8.088907739946926, 7.962488320214104]),
array([7.2674438325546555, 6.864833734600146]),
array([6.235509073802364, 4.777374754507105, 3.5773847129306824, 3.599843221358911, 3.2647053233957006, 2.598312514265663, 2.8607319857638496, 3.079074098924251, 2.4622729402041736, 2.3212404086490417, 2.1116484666577526, 1.8097315835850205, 0.9350997451446607, 1.740187058827558, 0.1782162488343969, 0.0]),
array([7.308160545627094, 7.45681545435209, 7.221142231848174]),
array([7.422730327804262, 7.446439913264518, 6.958218198812053, 5.322921672349543, 3.8303612251137653, 2.6258373664313917, 3.153468294858168, 2.816733140660012, 2.3499528293420555, 1.5180851511139433, 1.0647465856902798, 0.0]),
array([6.347104499725818]),
array([6.9724501525368465]),
array([7.218798991516282, 6.1051719665193405, 6.17434112950392, 6.758495252804773, 5.209267828148149, 4.592337082234725, 4.015510707359346, 3.8370228627450125, 3.4154169562223355, 3.4363968726407292, 2.1785271330144433, 2.14772262303782, 1.0333298705429672]),
array([5.125211849340454]),
array([6.6859765493532635, 6.683475424132594, 6.783353223709534, 5.54969442495288, 6.3028393243590575, 6.53733323357675, 6.491219398845474, 6.458781608130798]),
array([6.375208624076306, 6.453905038017314, 6.145245825797316, 6.593979028189899, 5.740895919132081]),
array([6.034671188035797, 6.326479000834415]),
array([5.97433356075121, 4.4198552488205, 3.6276163316438526, 3.1061663713113767, 3.109279740356986]),
array([5.940837230520108]),
array([6.48846111362595, 5.885501092207656, 5.542526015571033]),
array([5.771019749820106, 4.951620607225387, 3.6339471871752504, 5.228973005094172, 3.3858140930695932, 2.0760240918765938, 1.1133913900659769, 1.3156029185459186, 1.5452221249732567, 1.3932785120924, 0.4832429923016174, 0.25460256035803275, 0.7149369842777723]),
array([5.279469117996868]),
array([5.8205036010235816, 6.175972572968264]),
array([4.872879559336004]),
array([5.738812961060136, 5.66905819654577, 5.098651531323857, 5.316772796105449, 4.723489776682473, 4.802281820762035, 5.0193533235859835]),
array([5.408900387797663, 4.762370291169299, 4.577050776224669, 5.2576972179894135, 5.049090831403064, 4.297916508475659, 5.180186686065408, 4.753519509203385]),
array([4.270162429639465, 4.864280508827993, 4.653035601574243, 2.9822998839831083, 3.182248266521974, 2.9961992442058527, 3.168702749793911]),
array([4.179990522425008]),
array([4.5375150481190705, 5.0575723726404735, 4.112063577575552, 4.9235880857994045]),
array([4.539402017647964, 4.273851445059208, 4.065245148190079, 4.170066533282242, 2.703862320452177, 3.4583369466160168, 3.350228305750623, 3.358496661168248, 2.5824104240137657, 2.9631894910336323, 3.220057747239629, 2.9392795485646306, 2.740700287836631, 2.7623159996810376, 2.963161092205271, 2.7686880019216997, 2.266287010797728, 2.216997299768387, 0.875982850370186, 1.637959826539044, 1.0673394364847981, 1.0263453254209955, 1.5749841500400268, 1.295842004100943, 1.303836927616703, 1.0255800498715448, 0.9981892382009828, 1.597437962560756, 0.47680424958327455, 0.5389237915610415]),
array([4.287329016843136, 4.32394059741619, 3.2094653787682335, 3.5820654705710933, 3.258064726828046, 2.78041152401578, 3.507784079958113, 2.0638278006578816, 2.2907275408416155, 1.8883872902588306, 2.1724413312434723, 2.110713934252678, 1.493530377736973, 0.9264039587375247, 1.5929157366292224, 1.3100433524569364, 1.0919778033002485, 1.411029714986543, 0.0]),
array([4.526492387970473, 4.121041577028353, 4.0443950864421385, 4.448895143395522, 2.5812403237920245, 3.5591937083475718, 2.7645813248441127, 3.266382619449399, 0.8906437103200324, 1.6372973190436668]),
array([4.059784469675132, 4.1530130624833745, 2.9778022290758424, 3.4908232305699194, 3.214860084103632, 3.425364729079823, 2.62387702560833, 3.366089890399712, 3.460283940162177, 3.408434055499821, 3.2966040461618356, 2.8249264490916457, 2.9237367693396203, 3.2917994011343197, 1.9269452913033822, 1.9646291153762871, 1.9264262979582663, 1.8399148703038142, 2.3228457833899236, 2.3138730337054376, 2.045144195339942, 1.6981472097777017, 1.6992615776071391, 0.7989014407015669, 1.7951111074792259, 1.4201779840642528, 1.3337633616091504, 0.9974894519458689, 0.9885505975587592, 1.5408488033029246, 0.5916661364113215, 0.3537430464065209, 0.46352171774273465, 0.28214854673421197, 0.31589041858498584, 0.6251002374661092, 0.7246409468683098, 0.6607431364105449, 0.7303438234612909, 0.4521807320037385, 0.2729237191845224, 0.0]),
array([1.5859351997555626, 0.0]),
array([1.4044907282075905]),
array([3.687356466527114, 2.7986684496973195, 3.57565458622521, 2.032715572590256, 2.557134842844396, 1.5647942627248221, 0.9149297776228759, 0.8591849368064752, 0.8435648800576643, 1.1389591886220312, 1.1446219014275862, 1.2616197726452931, 0.0]),
array([2.6095685672988416, 2.6011935626547102, 2.8826119969426527, 2.0900550143080543, 1.8576745237703514, 1.8607635095587804, 2.3789390165154294, 2.549212791443705, 1.0534807380432396, 1.138921343573641, 1.3236947270793853, 1.4358360787895175, 1.5764690175000702, 1.1900210742387898, 0.5249594500660639, 0.36927196600666584, 0.6706229981150191, 0.0]),
array([2.824765352469856, 2.8160237574638423, 3.0826821165081584, 2.26453724533379, 1.8462536890212702]),
array([2.8624547286204427, 1.8462490580923157, 0.8780710578112821, 1.513506534270034, 1.7151452471157367, 1.0170861698219325, 0.6067279794556272, 0.2677157281723881, 0.22090774399080382, 0.48410455141594344, 0.36093132592338806, 0.0]),
array([2.0572352932365643, 1.629598478957185, 1.6367160435675556, 1.2337889510703874, 1.7423494557781118, 0.7888850350011729, 1.5177832773355227, 1.305379921608946, 0.7843244280031683, 0.9345642035888628, 1.4320905570459816, 1.1890041987696855, 1.2268875723325885, 0.7810238916156247, 0.19725414089366444, 0.45432956148206827, 0.18029262312183447, 0.7635724400389838, 0.4536667940207823, 0.051370835861758, 0.0]),
array([2.5976210268476168, 2.595177071693416, 2.164131750655662, 2.036282120225573, 0.8052538188150911, 1.254814063613602, 1.2475946736738004, 1.6604543155764446, 0.0]),
array([2.3808911413754266, 2.183563384959332, 2.2903508229758485, 2.1544272342218065, 2.1391783169127465, 2.427145440462979, 0.8063866602759306, 1.0911356302218365, 1.0763141422454305, 0.9011543843435638, 1.570654246806308, 1.0514602921096219, 1.2676865056683893, 1.6922413814531252, 0.7655209606303035, 0.6494319535213784, 0.3921572066907308, 0.5194861633550802, 0.21794543047605885, 0.5132282718490492, 0.0]),
array([1.5640995280671464, 0.8459703537894687, 0.856593942611055, 1.4980686093718987, 1.613334962164976, 1.1853196778221902, 1.4732807175602471, 1.0962903761117526, 0.7957478227863608, 0.9178424123328696, 0.658401671532936, 0.4318066721276686, 0.0]),
array([1.9368209703608306, 2.0073260999089446, 1.855861702793521, 1.5339911547982283, 1.59782864487906, 1.5848930450912442, 0.884045612245593, 1.0110834049077464, 1.581616097184697, 1.3066833763909265, 0.15391176137592177, 0.6335781670704662, 0.5026063976340014, 0.6782347969470638, 0.505588523140323, 0.40446120577785905, 0.6510414213696322]),
array([1.2723416380155927, 0.8394483115749065, 1.5024869273783232, 1.5808610497345066, 0.901225389871171, 0.9600762303856191, 0.5865704312594036, 0.6645139436077889, 0.1632934770410196, 0.14466750663204686, 0.1911387184606136, 0.37783083631473235, 0.391385724676015, 0.0]),
array([1.2800149701019539, 1.6182577581622335, 1.4077636462789174, 1.087109238388135, 0.3726630916726553, 0.26513415064434964, 0.49378816998475517, 0.4410396747118346]),
array([1.535610851363214, 0.13058877095052412, 0.5245403117887328, 0.0]),
array([1.2984609333795056, 1.2152626557966681, 1.5948050859343408, 1.1619014671724304, 1.049195976585656]),
array([0.8710327766342784, 1.0745328671461252, 0.2478203505540646, 0.18820184551671626, 0.5079061952833832, 0.5647611736038959, 0.0]),
array([1.0272599514362621, 1.12819773529345, 0.9798465963167927, 1.1644547423780265, 0.968529082013718, 0.6240163316906033, 0.16479460923852218, 0.0]),
array([0.9266261606302169, 0.9330984256260708, 1.0706390022527381, 0.5887589628203254]),
array([0.6322233041091324, 0.0]),
array([0.6147729238421669, 0.24694740253951364, 0.3024997575791012, 0.0648817855402136, 0.012009785449172403, 0.0]),
array([0.4957963850147683, 0.4952116211246871, 0.0]),
array([0.049993550582571364, 0.0])
]
d = [data_1]
names = ["39"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T20', 'T22', 'T23', 'T25', 'T26', 'T27', 'T29', 'T31', 'T32', 'T33', 'T35', 'T36', 'T37', 'T41', 'T43', 'T44', 'T46', 'T47', 'T48', 'T49', 'T51', 'T52', 'T53', 'T55', 'T57', 'T58', 'T60', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T84', 'T86', 'T87', 'T88', 'T89', 'T90', 'T91', 'T92', 'T95', 'T97', 'T98', 'T99', 'T100', 'T101', 'T102', 'T103', 'T104', 'T107', 'T108', 'T109', 'T110', 'T111', 'T112', 'T114', 'T116', 'T117', 'T118', 'T119', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T128', 'T129', 'T130', 'T131', 'T133', 'T136', 'T137', 'T138', 'T140', 'T141', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T158', 'T160', 'T161', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T168', 'T169', 'T171', 'T172', 'T174', 'T175', 'T177', 'T178', 'T182', 'T183', 'T184', 'T186', 'T188', 'T189', 'T191', 'T192', 'T193', 'T195', 'T197', 'T198', 'T199', 'T200', 'T202', 'T204', 'T205', 'T206', 'T207', 'T208', 'T209', 'T210', 'T211', 'T212', 'T213', 'T215', 'T216', 'T217', 'T219', 'T220', 'T221', 'T222']
def get_taxa_names(): return taxa_names