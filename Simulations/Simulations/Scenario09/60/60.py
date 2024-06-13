#!/usr/bin/env python
from numpy import *
data_1 = [
array([34.19625813322526, 33.59713281337515, 33.07435859165776]),
array([32.00693809950186]),
array([33.48124545405328]),
array([28.920042458861708, 26.906488859869967, 27.105361635090695]),
array([28.540870907028168, 30.199425672518956, 27.96471605915527, 26.937162736075347]),
array([28.559028202337622, 28.55994110115807, 25.695409150852306, 25.7928797000852, 26.294433825817578, 25.956880956086266, 27.033432551641408, 25.59088943770167]),
array([27.007243481256783, 24.760738389505146, 21.055952420211447, 19.433163408360848, 19.342696455031522, 20.393846333396137]),
array([24.21615312309602, 23.12905470642629, 25.28481968340302, 25.369295918238844, 23.52311733446869, 24.071109905611518, 20.699133842493605, 22.42374480052288, 22.949850289743484, 22.423181480002967, 19.645917507411056, 18.29049720154841, 19.706456272836505, 18.444308966576916, 19.880049897711864, 19.735686840926245, 19.60103190993683]),
array([23.783789559200667, 25.495519829792414, 22.682541450089236, 22.411427002207017, 21.2442946729905, 21.830143478370424, 22.882851606634528, 21.594999771328027, 21.367720286853064, 19.487408548922183, 18.59705636069856, 19.69326500119177, 20.00010592299693, 18.526314343764707, 19.344802832212753, 19.160454199828646, 19.293374224975004, 18.43220367583219, 18.481663865195046, 19.044234200316424, 18.64107010360262, 20.374225184924505, 18.266476963119143, 17.778435313282383, 19.316158012961836, 17.630078042041802, 18.842723392927702, 18.9045082568208]),
array([24.685012227317717]),
array([23.245473132539654, 23.530059901235486, 24.352911405062923, 23.694543551363328, 22.924604400217643, 21.86728651570626, 22.16087557756971, 22.725898470702212, 22.285890980195507]),
array([22.020203401587125, 17.349519906162776, 20.207029340622995, 18.017085473796573, 18.327623482693703, 20.37986658784909, 15.235033180091959, 15.394271476481885, 9.230886057189885, 10.566823721078743]),
array([21.55520797308821, 19.808446007009962, 17.051777386044016, 17.22836271364758, 17.262939991865153, 17.39300382093626, 18.875593394579624, 18.10953125189991]),
array([19.819002579213322, 19.30195088347904]),
array([23.23132992477816, 20.659315620450705, 21.18660139717065, 23.003705046980027, 22.80594438218839, 20.90866357795304, 21.377572389320974, 18.81386356004324, 20.158912207455508, 19.552534371518124, 18.97816088450931, 19.78961654601439, 20.161714379883513, 19.393404551497042]),
array([22.072089951798212, 22.388441583386243]),
array([21.706046408454885, 22.42998325678789, 22.464848833530983, 16.1319728283095, 16.459396150717406, 17.80810231047747, 17.36503550615549, 18.681491345868523, 20.150634321363256, 17.672972977529728, 18.298242177065763, 16.67868092213768, 18.020299890610644]),
array([21.53856587785097]),
array([20.542597457548602, 21.15013400817957, 21.080245460702965, 17.69074390334134, 18.808510487897422, 18.660037385933173, 17.142513299120548, 20.035485944666206, 17.223731732601433, 19.263684502778972, 17.37118666644673, 17.75324132907586, 19.200696820140664, 15.987177223599137, 14.097856975188588, 13.960303126557488, 14.527641195923199, 15.218203546618962, 12.33759421185246, 11.064042563922024, 11.396842915825927, 11.397355004363076, 11.44552793307127]),
array([20.75933084172708, 16.824662671138267, 13.85948476026466]),
array([21.622560045839254, 20.772203482360517, 20.163918877598913, 19.97310722433007, 19.851576290802583, 19.959406325659447]),
array([20.613501454105265, 20.938477450799336, 20.707914416797546, 20.483420581566808, 20.39988927802466]),
array([19.433122934996504, 19.65053907804454, 16.67598952273403, 18.381542352640523, 19.981549996552125, 16.63735399284775, 17.409468721680508, 17.8730754750854, 16.967198706441852, 17.904846431759232, 11.82398880885722, 11.846233255036642, 13.30675100719409, 8.201841521905493, 7.692189597538512, 6.444615043505336, 7.092458051345512, 4.218325474151564, 3.16251894799845, 2.360404010337744, 2.1461434979141463, 1.88153099919779, 2.2623139389138314, 0.4046205318260719, 0.0793231593655854, 0.0]),
array([18.590705025925658, 14.355988167358126]),
array([19.604932200467122, 19.934261280409693, 19.711376778405782]),
array([18.361763246786268, 17.00480770996291, 17.261354845062556, 17.815743041508245, 17.99018000402101, 18.354987091441412, 18.4523205510779, 17.471344391950566, 18.169722782027, 16.832631777862385, 16.04097713189325, 16.930960952269192, 17.5786346295909, 17.9694345274232, 18.258814000894276, 19.148880105332527, 18.13296750247133, 16.80455140942099, 15.3771769393944, 13.907439288319306, 15.011649478542852, 15.135542826194442, 14.219649252171621, 15.338212921826479, 12.911435963656215, 11.972300386003356, 12.08196217546222, 12.400117857476037, 12.018030133190068, 13.040530533369589, 13.529748789761426, 12.419742221366947, 13.219036024118365]),
array([16.749230648662245, 17.990656161843443, 16.96665743230627, 16.63943922004887, 17.86304829084727, 17.76884583661888, 17.35895881162046, 15.463855293973916, 15.25947827349132]),
array([18.44165254472164, 16.160754312906406, 18.21170623222424, 18.750717169525906, 17.767721676630355, 16.694933370764456, 18.929290052139628, 18.32141114163582, 18.907721444310273, 17.20649426248729, 18.34718267881947, 16.725844116649228, 19.06182045011395, 17.974502964399317, 17.264796321612483, 16.95912823629552, 15.43607965431266, 15.950195816764666, 15.705819864047026, 15.475492083032945, 15.729040279997752, 15.35570394838852]),
array([18.536639520481895, 16.83830674159971, 14.683689924510439, 13.957604653049158, 15.495074223897062]),
array([16.45751461281231]),
array([16.820824086201213, 17.43541543488376, 17.422356315231056, 16.175779670286975, 16.52822234085019, 17.082468266074507, 17.45655719299114, 16.110848703603278, 16.981707123857664, 16.54725321585954, 17.76947639070372, 16.14117999323479, 17.839024581507704, 18.132143546013495, 18.709228804764226, 17.976400677871222, 18.43653203861025, 16.7892070726179, 18.394438062750307, 16.86331196489554, 16.146498584234397, 17.506940515603436, 14.368058471352931, 14.976663634367084, 15.25364896123316, 14.647825335804407, 15.63810590230865, 15.767611731017576, 14.768983858858732, 14.63499628220559, 14.25639450978807, 15.96079215867003, 14.058192756737373, 13.06321677905115, 11.959658059736013, 13.107832483250302, 12.200902240157049, 13.69552897460699, 13.406608438026936]),
array([16.424828752820797, 17.810855283615382, 15.998111831466652, 15.938865071246036]),
array([16.63804202281499, 16.490326042765528, 17.326461165970986, 16.19412916457262, 16.15463310957868, 17.20138980969657, 16.249059132492953, 16.139908305666168, 17.68173154922436, 16.088033101635876, 16.335630711285823, 16.10794694034008, 17.70875744334386, 17.588609244477542, 16.487090814421247, 16.192032537570583, 16.72075042077772, 16.635852682569652, 16.583642768686882, 17.63091001051679, 16.22210656660217, 17.960054819923382, 16.626037123683655, 15.976604076364726, 18.26829515973003, 18.35579836917014, 16.193689263674994, 16.991183205847072, 16.313510004015153, 17.84330494134869, 18.167377202948206, 16.43887619362117, 18.152356439373335, 17.869398195165022, 17.155965527135013, 17.821732402043367, 14.210257098226077, 15.87961998306537, 14.802642386599567, 14.682422288152175, 15.313955546228406, 14.46299762919277, 14.404315324516231, 15.65945338351805, 14.00867759678315, 14.580704097057872, 13.852016798887927, 15.73390824259347, 15.162401400381299, 15.030834410796544, 15.293614269323857, 14.398999505258812, 14.042233777400341, 15.160301633515068, 15.3071625159265, 15.937573949324445, 13.85161694470776, 14.100458326917305, 13.874416326522113, 15.13211561387599, 15.740321547938752, 15.297393727685456, 13.527731381744124, 12.381117674342343, 12.309025935208451, 12.510709174433297, 13.049924752631563, 13.662824726125919, 12.848135932992408, 13.193193533469206, 13.526634373935627, 13.422887163828815, 13.429806144473677]),
array([17.066376802575025, 17.625751734235816]),
array([17.316549568904392, 16.903952848859976, 17.068981121218492, 16.98846137653962, 16.757005947303377, 15.852588483253655]),
array([15.973554168304801, 14.140299053783593]),
array([16.609601382282513, 14.595555865033033, 14.0757670557838, 15.696302130382254]),
array([16.34982323049777, 16.22295960907674, 16.433976077153428, 16.428362997387776, 17.263358386476906, 15.518620410924285, 15.288986512344058]),
array([16.188089967485173, 15.281175563560241, 14.66725878956439, 14.725895172977122, 14.441453722353085, 15.154830244659605, 15.053652813128373, 12.151501914304289, 11.862561024184158, 11.381252989380998, 10.299898848692148]),
array([16.193613106166435, 15.930650807069599]),
array([16.046892405239696, 14.092985787420815, 14.792076283646201, 15.637607187209424, 14.157969908042565, 14.853542782453154, 15.04087551822028, 12.044802152061685, 12.901577900922684, 13.096352653277544, 8.53028268349606, 10.963914497580468, 9.558515378827376, 8.557724902602862, 9.366492495314226, 10.815669431522887, 8.532037632875166, 10.019796747630272, 11.207913505653993]),
array([14.87338446111848, 15.214334921095377, 15.405067501942819]),
array([15.698489112295094, 14.850320507869418, 14.018866822577529, 14.120675339857268, 14.727351375113884, 14.490894847235518, 13.899893694877624, 15.149174581076315, 15.271888529806517, 13.99541591928239, 15.628209467813717, 15.343696048923304, 15.181137799737245, 14.682522568211175, 15.271171750855666, 13.992981724889828, 13.794873829926455]),
array([14.12243762647413]),
array([13.851936436631703, 14.539965140871303, 14.317802799776684, 15.162777245501841, 13.519646707815774, 13.017148744094863]),
array([14.97670605117395, 14.821419759843089]),
array([14.598450418733297, 12.98226003891437, 9.228545192705099, 10.081108700357376, 10.693279084957458, 10.462637614943205, 9.462562036844657]),
array([14.310897227712774, 14.609269587870127, 11.898740955215446, 12.06504100170572, 8.392563587337705, 1.8715926039208344, 0.0]),
array([14.538015780044564, 14.319150016454092, 13.12164409232663, 12.287368724891104, 13.512278889719138, 11.948023848262569, 13.260252260703261, 12.605070781567589, 13.248489688998863, 12.489620234726141, 12.064181910982269, 11.04300450053481, 10.511006685366134, 9.028898511272482, 9.750981870066754, 10.095658656140287, 10.286157846288354, 7.3090432671739585, 9.597854260414499, 7.443524378667871, 8.640313222047242, 7.736697355286645, 10.742360936739116, 10.677603953644056, 10.319550849236261, 7.65841159183282, 6.128487704786906, 5.336091210541992, 6.694717677903983, 6.743239364528065, 5.588512552895548, 5.8921166024910105, 6.727682783063622, 5.257468316277386, 2.2013768522390817, 2.375250995531087, 2.2976809096233115, 2.4682886896243716, 2.1396055459141166, 1.906301971976783, 2.210081529707508, 2.0347466335547892, 1.0730171269126383, 1.2099162594734438, 0.9122669067559668, 0.8037556665098611, 1.2277901094871444, 1.2996707055545902, 1.3512857038244412, 0.9870903229259556, 1.6459695667978442, 0.22838750130507246, 0.15509626136154797, 0.1546631249388004, 0.12170609234952036, 0.0]),
array([11.882694582378473, 10.166380151494154, 1.1918495809114629, 0.0]),
array([8.651544896501832, 11.276023322066363, 2.914927778377251]),
array([13.91702529203899, 13.978476384958663, 13.779700400785039, 13.303075426134052, 13.721591493753516, 13.673366951814494]),
array([12.303538880719666]),
array([13.381206459090834, 13.783940975251612, 12.678987965938354, 12.211522454142171, 12.57721837177743, 10.676990738708973, 10.848388456144722, 11.494814246779013, 10.849732014795729]),
array([12.022265750425737, 13.360010561526522, 11.260762686806707, 9.83583317548488, 8.813361749822635, 4.16404604259893]),
array([12.636367449746412, 11.685500061887565, 10.205316326236046, 9.661764852027073]),
array([12.574630425591927, 13.215820115653147, 12.903773716194475, 13.43042399614298, 13.005655877478029, 13.332621063333335, 13.31754832493794, 13.409738202429084]),
array([13.059847855204037, 11.680740724072322, 11.671131176382788, 12.733418606205387, 11.77065401562076, 13.077540412503943, 12.699830119425943, 10.574274167204338, 11.409631625710162, 8.098441651768276, 9.80675029535111, 11.544650080790852, 8.839990682254134, 10.280679023408652, 11.191377848580034, 6.689699957181816, 6.063992667718141, 7.049602559573161, 3.7141338399829786, 4.783222760243687]),
array([12.92247066327095, 7.791178430596668, 8.853703798294829, 9.116750568114366]),
array([11.87293988676309, 11.33797642649428, 7.677859758539845]),
array([10.058386306973317]),
array([12.23377607946061, 12.1199105479056, 12.081159871111376, 11.910627184090613, 11.413553831152543, 11.43256646089336, 11.022792917864793]),
array([11.715044345375517]),
array([10.489304145533753, 6.145613411401254, 3.584835941498552, 0.0]),
array([12.269218294195978, 12.084594949253413, 7.351231031649429, 7.542525321042693, 7.8403861764452305, 5.5353614434779805, 3.711784639718824, 5.145787721145228]),
array([12.162557618765307, 12.130761362484, 12.165228645390735, 11.425671873163125, 9.74178772448596, 9.517026576001808, 11.593539188301703, 11.494411243482796, 11.01639723742439]),
array([11.005705053718792, 11.173317081957977]),
array([11.835529933594945, 10.702651567780656, 9.064709089844637, 8.888944542735876, 8.111673124601513, 11.33343906974719, 10.221475194919602, 9.891464852256805]),
array([11.379503469316486]),
array([8.337466111575596, 8.823523120355608, 8.67056301701533]),
array([8.94653745307735, 11.053364621223423, 8.016629423266163, 7.395304469076252, 8.555402888406727, 10.500393651368555, 3.9241690153056243, 4.745180855439589, 1.977579546346343, 2.4366175093249485]),
array([10.690438946960574, 10.481011066755018, 10.85535500886083, 10.551839933065311]),
array([10.939161698309656]),
array([10.418730356407844, 10.545121487929924, 8.505369358717102, 8.163342972794382, 8.886859904602858, 9.705312096037508, 9.687750464597398, 10.227555431045417, 7.569457903870527, 7.737799852295215, 9.289778846856173, 6.231237762003385, 5.3260335386037285]),
array([9.82829920717311, 8.39426321830714, 3.614700118358085, 4.625025749172996, 2.4797479354506606, 1.3327810394433888, 1.106911163093623, 1.657739055147194, 0.02343716425448121, 0.049203312037727986, 0.0]),
array([9.523077583636129]),
array([9.630672887466156, 9.638403264042411, 9.39765220231166]),
array([9.170137580826843, 9.853596393642338, 10.61405881119091, 9.168958818102249]),
array([9.989213210451771, 8.105864042160029, 10.314238547239546, 9.975728391311767, 9.243261991254126, 10.13305030804852, 8.965224276735302, 9.899963227181104, 6.39656065382477, 6.276600610442181, 6.181199916453855, 7.007680772080886]),
array([9.930329142855385, 9.398827269964015, 10.161816087155504, 10.021066022001106]),
array([10.006883239014778, 8.278181732403603, 7.348125821591941, 7.405427255456489, 7.345812119992329, 8.174065570516781, 7.892250421971182, 4.5488505757326765, 4.482493385534861, 4.924513409409332, 5.074648592784566]),
array([9.857524546063194, 9.032204015818301, 8.922072551235479, 8.418624437921647, 9.339813207037514, 5.178422089291149, 3.9859343257482935, 1.6211856285430601, 1.0926685920519215, 0.5597752071109015, 0.0]),
array([8.696185636208352, 8.420097123134527, 7.829366160050824, 9.700418963493123, 9.48527901403077, 7.496860298373285, 6.475194474731583, 7.097089794928385, 6.9089711227018835, 7.14384065086163, 1.6788114686279798]),
array([9.759766759902941]),
array([7.952778183243425]),
array([7.768694678494215]),
array([8.148349112040403, 8.391160283402929, 8.599409167020573, 8.433238299690505, 8.219346102859735, 8.071377076844389, 7.828556692098305, 9.138123082913362, 7.583719514734362, 8.668943351693738, 9.21583209861632, 9.289980589486849, 7.587145150528111, 7.250261310909831, 8.88315928572919, 6.567774300298477, 6.674204582813816, 7.207528352796102, 6.66032958209569]),
array([9.062423540820548, 9.095200706785667, 7.764942839214452, 8.651519155195793, 5.906267597554406, 6.841929273152367, 4.061913883526586, 2.1715747379039514, 2.0986577460005047, 2.5777576338446373, 1.335157721097558, 1.7230061242227621, 0.861217926923758, 0.4910373896723086, 0.0]),
array([8.683650068145736, 8.345223215028433, 8.093948048603261, 8.431050408051062]),
array([3.1650411081970713]),
array([8.18371463956396, 8.301547756614093, 6.331786991189522, 4.858104696199957]),
array([8.370697968669495, 6.4401318260629115]),
array([8.359154911372483, 6.347739938109207, 3.6283139606622465, 4.958944964145903, 3.1519145873462073, 1.9556061671408687, 0.7949598169009597, 0.13898795093431449, 0.021418328076990192, 0.0]),
array([5.681576805975356]),
array([6.3120470551049035]),
array([2.60335828777178, 0.0]),
array([5.097860651505303]),
array([6.713645250158813]),
array([4.447569481714957, 4.869089648588979, 5.2173781187985036, 4.976610782394898, 4.748062567538609]),
array([6.302549752405421, 1.7903348109011652]),
array([5.787431527376208, 6.500857580631187, 4.627062631489545, 5.094483661126962, 4.900648408723788, 2.1149353531385144, 2.43222623186391, 2.286916489550954, 2.0123770150737768, 2.178785815038973, 2.0927359353083017, 2.4220327247629356, 2.413002871978403, 1.9608127879101278, 0.797715117683232, 1.213365302065296, 1.417401405249437, 1.4913476324352877, 1.2367655310238432, 1.0439540298006502, 1.206477689473195, 0.9875998870569047, 1.333341160933457, 1.1643774097696937, 0.3202989875352251, 0.7636516315141584, 0.7079900304546719, 0.2558722816552178, 0.1022161367335348, 0.0]),
array([6.294960865335026, 6.276040801706839]),
array([5.6208620339326965, 4.212708891416509]),
array([3.2894171853829923]),
array([6.30833874058946]),
array([3.054064831636707, 0.0]),
array([4.871236995659638, 3.1926237576088137]),
array([6.249854150392126, 5.788755607503389, 5.2674148076501055, 5.217160992284954]),
array([3.948783458756413]),
array([6.154213190299446]),
array([5.394232611622984, 4.5724577341439465, 4.890726584981787, 3.859807678365685, 3.930195970696049, 4.280758286990597, 2.8881557890253067, 2.4100732781487975, 2.363258697373878, 2.339389496660724, 2.546955305613694, 1.8233552179828707, 2.5587083333849017, 1.9624856538924775, 2.125193331228422, 1.9565397292104023, 0.9467588661065655, 1.420821289640892, 1.2898919933182484, 0.9025205147424976, 1.6825749091826667, 1.3517812892082852, 1.3128980552584175, 1.774615992461415, 1.4834288138146512, 0.7135783041198719, 0.2423400471456787, 0.41744558355803946, 0.03651920571746274, 0.04856197125051946, 0.11803171544628537, 0.0]),
array([5.293943371315057]),
array([3.695990415430437, 1.6295525239985449, 1.043992312005337, 0.12297335285680158, 0.0]),
array([2.0073502206941822, 1.7743415988431892, 1.0702553857456758, 0.0]),
array([4.179723509131446]),
array([3.6109591767234477]),
array([4.194439480844061, 4.465247932590341, 4.736851584007916, 2.802491324984053, 2.46819231633903, 1.9657196217354829, 2.430200062876505, 1.0044363676629258, 1.0943500054464914, 1.1615624910758782, 1.4627888229230526, 1.6457511413183095, 0.9587357564568131, 1.4882970752691829, 1.4936275866993487, 1.7421241178310671, 1.1633885041634162, 0.15028291959649864, 0.056206124544852284, 0.0]),
array([2.4750317667826907, 1.0392738524545448, 1.7786950416639309, 0.0]),
array([1.9471211427439519, 0.9808717817553102, 0.5476906974651627]),
array([4.29448136874291, 1.7471339109881916, 0.9635067276099726, 0.47348621807706387, 0.0]),
array([4.1531599723915695, 3.9633699354386294, 3.5505253660668403, 2.4879194031066145, 2.38744565065599, 2.3444594816808757, 2.4056375967523382]),
array([3.9864511346541125]),
array([3.569317118327204]),
array([4.096907825397596, 4.046416782834148, 2.1330857971146666, 2.5512121343028076, 1.7682637788322297, 1.4164578747213419, 1.3733725860168433, 0.9922676835798973, 1.395750721028685, 0.25053344673940914, 0.5230157248634784, 0.0668885114951667, 0.11716304210719214, 0.0]),
array([1.4731187543417172, 1.3733577903665832, 1.4591620021389151, 0.9495137698892693, 0.0653784708149042, 0.0]),
array([3.8976026168486753, 2.589685797169774, 3.2389275940444087, 2.7915730046480864, 2.579272285406894, 2.1711759915251907, 2.4970259664958263]),
array([1.8921705406703349, 2.1297146553419233, 2.041800682373829, 2.562716206770767, 0.9629552039746959, 0.810445415523081, 0.8582790309709586, 1.3464015105816942, 1.6081128899249393, 0.5636333493473307, 0.20559011373643288, 0.08235540915468711, 0.09422712930630826, 0.013973544366151336, 0.0]),
array([2.0052490286752436, 1.399263584703506, 1.2481204975970759, 1.6647834342804644, 1.1129259544937553, 1.5143978987064888, 1.218795938052577]),
array([3.031781815574921, 2.6449254125876425, 3.30964370805633, 2.7802595217101365, 3.289804712947534, 3.128219451322686, 1.9570037204446842, 2.5049561547408787, 2.3524158311992567, 2.5325543115715066, 2.130857733925676, 1.8414545859472962, 2.0263491659592994, 1.8251860972671938, 2.5797454887540123, 2.249858198754317, 2.2054964390572755, 2.2613418070876796, 2.5119386715682426, 1.8354617785300893, 2.408672655298462, 2.576313364769752, 2.4277756901395176, 1.163518842584202, 1.2821110856953153, 1.1966743327793548, 1.5434083099405997, 0.9128385437730339, 1.000281474344661, 1.0666424922979219, 1.1605349593485488, 1.6933941351411788, 1.3217397169000984, 1.0823212868553853, 1.54501322522041, 1.0022702661345648, 1.229489309788198, 0.9650249584922509, 1.3391797509281855, 1.4563574202738243, 1.7302515668524268]),
array([2.748133540433604, 2.430486172617063, 2.367886864077536, 0.0771446607686544, 0.0]),
array([1.2364242740982014, 1.0366254686427778, 1.5897552872638334, 0.06199073030252221, 0.09963660731878315, 0.054765514224635034, 0.0]),
array([1.967229774923851, 0.0]),
array([2.6231173609681364, 2.863983012195283, 2.372861651288835, 2.078556090076885, 2.277969211252112, 1.8023758973744277, 2.112115100013701, 2.145996836692895, 2.4990818622612725, 0.9877378645146506, 0.9498566303065687, 1.1448440297758897, 1.3556633975906442, 1.3782711291355232, 1.008418025377455, 1.7383364498533982, 1.37974874498707, 1.0516120903634687, 0.8169936241470366, 0.9787314144565836, 1.318840651788037]),
array([1.9597423073226627, 1.9102818056965818, 2.5768512471147536, 1.506403115912441, 1.6464661266885565, 1.4676908844893741, 0.5302674466925597, 0.0]),
array([1.8034056627796602, 0.0]),
array([2.196927178970355, 1.687372697272062, 0.0]),
array([2.8951504280054627, 2.4929108649094136, 2.204504717680124, 2.113588145443678, 2.0408611956895655, 2.5632147555953275, 2.5094828686917916, 1.8223947753954146, 2.4517317231078346, 1.9474242838834304, 0.9098007439356811, 1.3082648426824763, 1.1023430650828383, 1.715216947160818, 0.8463769758825368, 1.2367976089468273, 0.8816532505428386, 1.6881447279002637, 0.8957413159072822, 1.1281332699853808, 1.1897031219422833, 1.2115958713907466, 0.7911340029462102, 1.4644710364691231, 1.7225012566407396, 1.221323785973059, 0.8133877955930278]),
array([1.8182995159385191, 1.2632438205192953, 1.7744899072328781, 1.02635850004029, 1.4784806138144304, 0.9225959541134411, 1.6795264244804402, 0.7992865704754484, 1.1060790732669703, 1.4870821391657403, 0.8254585971607704, 0.2881323526977024, 0.2320466302950619, 0.0]),
array([2.233461315393752, 0.9023985301214618, 0.7524725828486308, 0.0]),
array([1.9929617432110802, 2.066453279243635, 1.9236072905886945, 0.8461208918136353, 1.1747920053939145, 1.772422618447577, 1.3298726132778431, 1.6582974104806767, 1.6378236062163742, 1.386128860607428, 1.6410238803480075, 1.2881781897400717, 0.9010052972516227, 0.609644188508191, 0.1343194249830747]),
array([2.1942526957453157, 2.0433529736512477, 1.3710900372517054, 0.9958648267848748, 1.448346341619846, 0.7889694539547172, 1.4984494017639318, 0.08190414863675669, 0.0]),
array([2.1441584560586913, 0.9838209347417741, 0.8581509598140804, 0.8760696853226212, 0.8431132265919034]),
array([1.2073574878121627, 1.4539602475646942, 0.8316981511098678, 1.1170816844298357, 0.0]),
array([1.609178400028023, 1.3194911639286673, 1.205097070787581, 0.9727686559821986, 0.0]),
array([1.895404423712368, 0.857568768828267, 1.6984589427659822, 1.6632788781289567, 1.6849932550597375, 1.1776655793415887, 1.1886479975159059, 0.04314743517140336, 0.0]),
array([1.9449471979969457, 1.539481142093229, 1.660953001646864, 1.1980552386777958, 1.27460673981345, 1.4729870811813526, 1.331344629172985, 1.723413304401818, 0.9181162186288625, 1.521346403015826, 1.3927151328418172]),
array([0.05132266300973305, 0.0]),
array([1.669007709612805, 1.7171170754205034, 1.491076313622036, 1.3576077377738076, 1.2463648073336198, 1.4454973996012777, 1.684679197929842, 1.3466332436318527, 1.3645933622315207, 1.7025358951627627, 0.8610311982673357, 1.6872268058292406, 1.7395916753889322, 1.4064723893452888, 0.8740543938678438, 0.07832909086806225]),
array([0.8485322502803605, 0.9830440591650388, 1.329872387827792, 1.39235915458505, 0.7192155845689833, 0.0]),
array([1.3397632348143405]),
array([1.515110995794639, 1.0443411965580398, 0.0]),
array([1.377571729018076, 1.222033045422828, 1.0412809223488875, 1.1783452776154515]),
array([1.3507368978268424, 1.0345247034150655, 1.2316976608907464]),
array([1.4763712426010667]),
array([1.1967537717157746, 0.0]),
array([1.197016852759765, 0.0]),
array([1.0804279199386952, 0.981186309958626, 1.1916684900764463, 1.0921211624386453, 0.9650316835179416, 1.1533650107899163, 0.9888325373616074, 1.1084797173496834, 1.2762290255009001]),
array([1.0063223500021206, 0.11178850873361788, 0.04363390692519245, 0.08285772202448462, 0.12215845797574999, 0.0]),
array([0.9766527488558503, 0.8375304892590263, 1.0435952965875548, 0.19776153524888018, 0.5330091235765827, 0.49790058846713514, 0.0]),
array([0.8863692637414713, 1.1340536014451421, 0.03577233482021909, 0.0]),
array([1.0824612637053272, 0.8979274661167711, 0.027544329508974147, 0.0]),
array([0.06087695088898204, 0.0]),
array([0.9155936752616549, 0.8146314321410939, 0.7571858455763721, 0.4019740581577587, 0.0]),
array([0.22931256238949016, 0.617209825604596, 0.6873206547265764, 0.47916303043041986, 0.0]),
array([0.7479491085239701]),
array([0.16532623223968435, 0.39000288611071454, 0.06950727532886951, 0.0]),
array([0.44230235292812603, 0.6082335821694059, 0.24114832052788038, 0.10950577203251297, 0.09077733630252693, 0.03769328046177624, 0.06862259480118696, 0.0]),
array([0.010516779586539987, 0.0]),
array([0.07791083920429251, 0.0]),
array([0.27394722306904895, 0.11741259346194037, 0.07471226257431182, 0.0]),
array([0.34815921379282005, 0.10770090558297518, 0.01889735530928012, 0.0]),
array([0.040466238826182915, 0.022110243810968058, 0.0]),
array([0.11412830553566966, 0.0])
]
d = [data_1]
names = ["60"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T5', 'T8', 'T10', 'T11', 'T14', 'T15', 'T17', 'T20', 'T21', 'T23', 'T25', 'T26', 'T27', 'T29', 'T30', 'T31', 'T34', 'T35', 'T39', 'T44', 'T45', 'T46', 'T47', 'T48', 'T49', 'T51', 'T52', 'T53', 'T54', 'T55', 'T56', 'T57', 'T58', 'T59', 'T60', 'T62', 'T64', 'T65', 'T67', 'T68', 'T69', 'T71', 'T72', 'T73', 'T78', 'T79', 'T80', 'T82', 'T84', 'T85', 'T86', 'T87', 'T90', 'T91', 'T93', 'T96', 'T97', 'T103', 'T104', 'T105', 'T106', 'T109', 'T110', 'T111', 'T113', 'T114', 'T115', 'T116', 'T117', 'T118', 'T121', 'T123', 'T124', 'T125', 'T126', 'T128', 'T129', 'T132', 'T134', 'T136', 'T137', 'T138', 'T139', 'T142', 'T143', 'T144', 'T145', 'T148', 'T149', 'T150', 'T155', 'T156', 'T157', 'T160', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T168', 'T170', 'T171', 'T172', 'T173', 'T178', 'T179', 'T180', 'T184', 'T186', 'T187', 'T189', 'T191', 'T192', 'T201', 'T202', 'T204', 'T206', 'T207', 'T208', 'T211', 'T212', 'T214', 'T215', 'T217', 'T220', 'T221', 'T222', 'T223', 'T225', 'T227', 'T228', 'T229', 'T232', 'T233', 'T237', 'T238', 'T241', 'T242', 'T243', 'T244', 'T247', 'T248', 'T249', 'T251', 'T252', 'T253', 'T254', 'T255', 'T256', 'T257', 'T258', 'T259', 'T260', 'T263', 'T267', 'T268', 'T269', 'T271', 'T272', 'T273', 'T274', 'T277', 'T278', 'T281', 'T285', 'T286', 'T289', 'T294', 'T295']
def get_taxa_names(): return taxa_names