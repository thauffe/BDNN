#!/usr/bin/env python
from numpy import *
data_1 = [
array([31.2603574821041, 30.352039457601002, 32.49196412190288, 29.838529091376305]),
array([28.183432543997867, 33.131559325370645, 30.689008732387816, 30.389119625304804, 33.72474781288111, 32.96808444410283, 27.87255835537073]),
array([31.261600385428412]),
array([31.49065823934093, 30.013703339003094]),
array([28.128586965760466, 27.18501937517394]),
array([27.403134535096793]),
array([29.054959073590943, 25.522927451629506, 23.05922929986207, 25.748856694096446, 23.70940956811831, 27.690737055740367, 24.55571947578122, 27.300423335207995, 23.442836056891643, 26.06770391260358, 26.014937817690157, 27.616489673627473, 27.712866419231958, 21.428246915514457, 22.852832212393114, 22.540134267098292, 21.19795134966614, 22.73723538325061, 21.814619943026713, 22.062957753020918]),
array([30.002820969382345]),
array([29.156287139391672, 24.997486978935015, 25.36370274126037, 26.33296038333958, 23.198999088883085, 25.857842237853802, 22.750375665553864, 21.216153868440767, 20.886376934087366, 20.722242697034915, 19.51292641361203, 19.954792411690732, 19.68811103178257, 20.33369447972613, 14.197343683319518, 13.675322033776286, 13.456264925171295]),
array([28.47983100707251]),
array([24.62582136310741, 26.896885012090056, 25.54767628730732, 23.343139588345018, 26.0362727772522, 25.361407692433716, 22.440036742909083, 22.01066391688176, 22.80432011681186]),
array([27.511404852315522, 27.458536482828443, 27.322898905001413, 27.5931188060209, 26.570494588627685, 27.093415234624068]),
array([25.836506968934543]),
array([25.864395823793657, 25.39582273711991]),
array([25.297549358773868, 26.173914378657837, 25.94283864013804, 26.29632175016021, 26.979995719368613, 24.28060369745587, 25.052309042484385, 26.232494576050538, 26.84669479708103, 20.70186986105905, 22.771178980214202, 20.68344557342307, 22.88509787007774, 19.969660470811238]),
array([25.59309675503694, 24.418439189699164, 22.012807391788694, 22.040967107067306, 21.69399455385871, 21.815821169166927, 18.20048165207419, 18.531109045236374]),
array([26.36192777507273, 26.278845371712656]),
array([24.985398458447186]),
array([24.689656084449364, 24.744225742517372]),
array([25.369166122335645, 25.38285688179861]),
array([23.711679852779113, 23.592781652802174, 24.52302470507724, 23.53503320435959, 24.163621051003958, 20.93661070527402, 21.05330155722673, 21.2401558768388, 16.246449509782032]),
array([24.110242968269702, 23.188338890137548, 23.90657528562943, 23.707402301692827, 24.55136686238208, 24.31650969509881, 24.069753635640225, 20.68733687351533, 22.92621099277065, 21.520208185574607, 21.867435670696203, 20.526739564368622, 22.70983422926039, 22.464160484889657, 21.800942027399568, 22.818881563411978, 20.086866722678337, 19.669225293511484]),
array([23.683924214170673, 23.480839881049043, 23.2949343416038, 23.53277524828663, 23.6862011043016, 22.640807244113965, 22.692640134538888]),
array([23.845658088088577, 23.94311505979044]),
array([23.655713625766488]),
array([23.455136020435727, 21.057870161538293, 21.65504867807815]),
array([22.945668005073315]),
array([21.944394759719675]),
array([21.02372428556606, 20.597549518834875, 20.91577102250024, 21.254342958185212, 22.547273185666818, 20.65709903517831, 22.66016929808209, 20.719981928017596, 21.699269871975552]),
array([22.079993998281978]),
array([21.82557219834871, 20.891424995847157]),
array([21.334812927491562, 21.688478444414272, 21.74273688121838, 17.664236674005718, 19.79045860982731, 18.710357862929577, 18.8906656920779, 19.989934390083505, 19.87399794324316, 14.44101398946528, 13.416462221543563, 12.93351436265875, 13.675012747689438]),
array([20.55155917941972, 20.641528055694693, 17.481316613654933, 18.133992175503145, 15.19284526764337, 14.604313266142565]),
array([20.98440877788173, 21.277994716251477, 20.96907814215199, 20.8019396912055, 19.506469154637777, 20.388816950156336, 18.44126977216676, 16.653792215351928, 16.12867853383207, 18.03496839594293, 16.57094767630497, 16.639108745138206, 14.178526173928521, 14.845084793983405, 15.395270036461888, 13.918045858413578]),
array([20.8502529610584, 19.869435289050617, 18.220571767108527, 19.16090500868523]),
array([19.60727680284842, 19.13683668422186]),
array([20.58908315318433, 19.10727254070783, 16.344207047819943, 17.66882877237667, 16.926685693981312, 16.995482054639872, 19.414309598698804, 19.632927735642788]),
array([19.828560021337783, 19.512558622549403, 20.13188228415467]),
array([19.63599224710215, 18.81663842253114, 16.461442635528474, 20.013900692557787, 15.569173202359284, 15.538228388515236, 15.266834024529382, 15.895863353661335, 14.512291250780631, 12.86647159921074, 13.034369131827427, 9.926126405266313, 9.962351555857534, 11.43540192678657, 10.15391510869429, 8.702906409100986, 9.08182186140999, 11.276506620175754, 7.856556241832674, 10.322710505877946, 8.406838247165185, 9.121681574256998, 7.11463497011934]),
array([20.46346834812703]),
array([17.415953191303373, 17.588471752079844, 16.91740619029009, 11.673614319264605, 12.81846924595915, 12.807783953093232, 9.537757611377044, 8.779824064091759, 8.873559425904192, 8.875191899265948, 7.8007270510493605, 11.42582055189104, 9.357328905177996, 11.136585739691018, 7.970978914559918]),
array([20.385716928042672, 18.77040130656476, 16.95405902405692, 17.805972366946868]),
array([18.421228174611738, 16.289696721654238, 20.027919361666584]),
array([20.075263331946715, 20.12919034572473]),
array([19.38336852182756, 19.393612641429108]),
array([17.064236352562666, 17.016414924277587, 17.43237789336139]),
array([18.940342695489036, 17.6625680215682, 17.887449597643737, 18.740976344041847, 17.613597245902056, 16.85728249842889, 17.784948077426485, 16.221270657519142, 16.303204361549426, 17.240905694207903, 17.863288045760648, 18.728732899914053, 15.56553482424004, 14.183475437574627, 15.685547073871867, 11.698961736178763, 12.963138715712915, 13.492237745885603, 13.37469230532402, 11.984159626537172, 10.676282482959245, 11.16693390041802]),
array([16.922843478835222, 16.728402063645806, 14.212814644527167, 14.731485540076129]),
array([18.96432218454399]),
array([17.76078245670428]),
array([16.535482936001156, 16.732159370392278, 16.373692274495617, 17.485936597911362]),
array([15.070066212289785, 15.690142572063158, 14.823300951703834, 15.016289548492836, 12.037354036323414, 13.185723583747121, 12.632298789222693, 11.06455626166717, 11.432736188854106]),
array([17.802055009277804, 16.025182431426188, 12.47789463007698]),
array([16.222432578449588, 16.260582642120493, 14.133320704264277, 15.002910631020805, 12.276489458071046, 11.731244933559474, 13.276451317410359, 11.649053367852849, 11.555181442884615]),
array([17.813768111219947, 18.571073928543196]),
array([17.515902116434813]),
array([17.973835906621794, 13.587405092671196, 12.429696137304054, 11.491783092373954, 11.503878451700011, 9.886643451143701]),
array([17.825354072824652, 16.893794388634834]),
array([16.626343076552562, 17.20735233211733]),
array([16.873769000714738, 15.836045710079922, 13.770308362232658, 11.433897807746046]),
array([16.96508835416644, 16.706537937125187, 14.437558876195375]),
array([14.680404559595464, 13.372638590753299, 13.81092674651445]),
array([14.803117928828133, 15.538259865761267, 12.97857405066481, 12.309575855037993, 13.735783190093397, 13.492833126639624, 13.227684909829412, 13.445612775446612, 12.817473130779375, 13.640355204723196, 11.396452386465873]),
array([14.902880330697883, 10.147404960559983]),
array([12.285601085843826, 12.068378687430581, 10.384064652849451, 11.37280721118274, 10.969711569273976, 11.499838481957042, 10.57695520171577]),
array([14.085164633110642, 13.663123584638768]),
array([14.279597578198906, 12.388164141237853, 9.238918970897357, 7.668472458086553, 7.290577600121194, 8.956840844650266, 7.454917049206724, 8.614208803070142, 9.691848007630131, 11.194436753069906, 11.243899935370075, 6.888673186301156, 6.726168627617282]),
array([12.533451260704155, 12.903832621256234, 13.678987105423571, 12.052724050544983, 11.745508042012519, 11.572102429249517, 11.471282419789663, 10.002705562790307, 10.87859210177563]),
array([9.338025263421928, 10.741891880832974]),
array([13.154341028522271]),
array([12.303026937515, 12.880666744943715, 9.357688868474177, 8.789306425028254, 9.197825542189834, 10.658342177104625, 10.516594843342025, 8.985097275195027, 9.208712703550402, 9.138325393588797, 7.761210143440763, 6.907643473135391, 6.112812167687744, 4.253341591243383, 3.6282410643407856, 1.4971855234925353, 1.6080941848697212, 1.6788290877258871, 1.1266143167008686, 0.9308158732331786, 0.9268788682497784, 0.34091206412816605, 0.201882558888269, 0.36255321336344176, 0.7275853277079046, 0.45768704778395425, 0.5076346557548013, 0.0]),
array([12.491613698159776]),
array([10.94492671968017, 10.985396248774151, 10.657821722027492]),
array([10.545025808742434, 10.737523913620475, 10.392380031113456]),
array([11.904574483315292, 12.465359468358326, 11.61071764564558, 10.505567956481388, 10.028344106573247, 11.464784357902795, 9.943511161849429, 11.533061553580344]),
array([7.472436529776114, 9.499729208295545, 10.808179564494392, 9.993466096955105, 1.6737613244657128, 1.7745407752592162, 1.7708105603215285]),
array([10.120275868140155, 9.669387071871759, 10.552450157540859, 10.384558896779065, 9.090041445586783, 11.233026063896922, 9.790033644087687]),
array([11.572403327892166, 9.058808636214403, 9.264859418858034, 10.749344922553261, 10.299918077602417]),
array([8.952137613270592, 8.58380560075243, 8.841900757563096, 10.655180001875861, 9.922370861050577, 11.13155207680283, 11.389817053173811, 10.104181104674531, 8.779283350014213, 8.771995997089583, 5.377472231901889, 6.456297274661127, 4.4933319528418165, 4.294572993026403, 4.404557658396909, 4.403506344746731, 3.407279182482031, 2.1788834295074566, 1.6997145784790435, 1.073173701687613, 1.499167515825141, 1.0472395281025624, 1.30073118419798, 1.1326034909518659, 0.8294692249269501, 0.9032905625461682, 1.6422281626977582, 0.610672473637791, 0.4884854438991915, 0.0]),
array([10.001119930240574, 10.578071124334057, 7.936583265425818, 7.491313714615817, 9.66355081409663, 6.2710089120550245]),
array([8.216283840734587, 10.23980773014937, 10.472299533178015, 10.633555139711667, 10.749517011692678, 9.836018985991513, 10.93078730043023, 7.564096466106567, 10.949074342636782, 10.928357582414966, 10.20994226496981, 8.811377108206548, 7.2082152281262895, 5.928407736312453, 6.80299954113995, 4.36619757415949]),
array([9.509630848212499]),
array([10.14042884969269, 9.767861087126443, 9.74301139271296]),
array([9.267855585227501, 9.490156847599685, 10.531973139426292, 10.535392716343503]),
array([9.959241236529934, 9.808382633622472, 9.531506243960747, 10.42113752435258]),
array([8.476914930712587, 7.93380088963148, 7.842814167973758, 3.657026681820167, 4.276123496638872, 2.0724404470599755]),
array([7.57482593227458, 9.309757463562784, 7.912966261531418, 9.806714808476633, 7.69136150499417, 7.306130579648581, 9.342601237897421, 7.198178486600823, 5.994254037802751, 4.568066964410029, 1.9142278479420956, 2.2723548715380297, 1.2889097298479246, 1.796225602400629, 1.140897776719624, 0.6920059846912245]),
array([7.975028993068404, 9.827842766422531]),
array([8.215716222985337, 8.30167668651971, 0.13315348686536232, 0.0]),
array([8.271113802450193, 8.74474404872096, 8.634038998540188, 5.689839439255664, 6.5045060078922, 5.947007019823319, 3.790781453114567, 4.597427031125413, 4.065159561440196, 3.0900946615405513, 2.6816867143185306, 3.4565541745158925, 3.0551239671318067, 0.896310120807563, 1.0078161479213705, 1.2960572819363922, 1.7447371258553206, 1.0821651021720122, 1.7077868850419846, 0.9342758950486447, 1.0483914924898021, 1.6423598988860912, 1.4654176416376898, 0.20712189900284816, 0.2271347338996682, 0.43657249866369136, 0.5754737146848264, 0.0]),
array([7.54501583991461, 9.612763064128504, 9.285311815930548, 9.93607666419145, 7.277152478214898, 8.314390239456062, 7.755605349952287, 9.838916494055535, 7.065011096505556, 5.558555402640404, 0.8583992927472829, 1.7830117175925153, 1.490297655286278, 1.706258887399793, 1.2662017427292862, 0.7833573377079392, 1.4246469061381621, 1.1852958111644976, 1.701243464997007, 0.7516630624343864, 0.6611475057612309, 0.272012732415582, 0.6684211726985089, 0.15914920348606787, 0.4315181347319562, 0.706419415643104, 0.0]),
array([9.342504373868056, 7.7816003454696245, 8.987682407910125, 2.7393650040693314, 0.8947923093397078, 1.225041257450381, 1.161890109728354, 0.39705166193072344, 0.4152006726204993, 0.45899641271047253, 0.0]),
array([7.9969643007410856]),
array([9.030433184599147, 9.258914080703796, 7.345933884740207, 7.816275305566738, 7.784899233205319]),
array([9.28546713627855, 2.960877274439107, 3.042515983607804, 2.9981253174575215, 2.011484570765077, 1.469356514220415, 0.30848951063157226, 0.6684506879052687, 0.24102780800349566, 0.0]),
array([9.453156810804773, 9.25107654104857, 8.287480447196108, 9.13697924962815, 7.5696672645637895, 8.549619517561592, 7.234970986150215]),
array([9.441737262093145, 9.131522123963462]),
array([8.454358139184336]),
array([8.085157455532531, 6.912043047565038]),
array([8.374798511737868, 6.457995810532084, 5.145237593237225, 5.004658064776288]),
array([7.940868986894445, 6.968670932404126, 2.9332377849044042]),
array([5.995575200439454, 4.327777940911951, 0.954840069662282, 1.3043582583846367, 0.6623184163778824]),
array([6.521844296969851, 6.228489024832337]),
array([4.068868960323827, 5.283144516210932, 1.3265440860029276, 0.7308012357640737, 0.0]),
array([6.2245503553853085, 5.088647121250415]),
array([4.242286488305056, 4.197513972979612, 2.9625336384376304, 2.4662524472452443]),
array([6.019080692389096]),
array([1.034286195654603, 0.0]),
array([5.5754133186694705, 5.964661377151096, 5.188694935437661, 3.3341338812838632, 3.3135175863685937, 2.197042447068454, 2.237557462744431, 1.2351720481114654, 1.207383280968724, 0.9269793182734928, 1.1167610548425857, 0.31038674896599583, 0.4470732345582557, 0.2549606285964179, 0.0]),
array([5.447924896337985]),
array([5.603696121137554, 5.693901052148597, 5.597431717676824, 5.856822212302605, 5.720209255251235, 5.127766737895498]),
array([4.468790604600891, 1.9130722325941318, 0.9612541675910559, 0.4120113104623542, 0.5424198131010589, 0.2349316383193769, 0.0]),
array([4.72140743273744, 3.859416159411179, 1.8562454008102784, 0.8681157740795301, 1.0087423669723181, 0.8420572283936508, 1.5981761406348913, 1.551941619889857, 0.6669111706124727]),
array([4.87519810163056, 1.980108357652079, 1.523161485825639, 1.3413704148784698, 0.8909675974478205, 1.534070414996639, 0.9594420009389236, 0.6374879763006065, 0.0]),
array([5.472439616800138, 5.454430001398541, 4.3154235655926945, 4.472509346175325, 1.3939196320049907, 0.6014744078832948, 0.0]),
array([5.324789116849291, 1.9029629082586146, 0.3929898836569409, 0.0]),
array([4.548212978440085, 5.036002317699656]),
array([2.9815738547788664, 2.8375377991553727, 2.070918407890958, 2.174659980606186, 0.8459929091778761, 1.2850422959766434, 1.5015004111650954, 1.0034735638967942, 1.4756297907569498, 1.5309245093736956, 1.0313639069361176, 1.7693569787191932, 1.3224539097780035, 0.554340072390983, 0.398191256074785, 0.2648018698648199, 0.0]),
array([1.2025633384199144, 0.45568203742935776, 0.0]),
array([3.627445228523837, 3.963372094127635, 3.585494833003646]),
array([1.0007305739952912, 1.674223456982127, 1.3236281684739888, 0.3439220586745536, 0.7309936202006173, 0.0]),
array([3.5818516490214845, 1.5660535841745522, 0.8395435000984552, 1.2457801347914637, 1.333544318016249, 0.5047630588736925, 0.0]),
array([2.6952107414617803, 1.9363561977080739, 1.674201202183273, 0.41097233561275953, 0.34077093267028397, 0.0]),
array([4.7796479945133]),
array([2.882227821608135, 1.1212273643396644, 0.0]),
array([4.226222034760377, 2.2522303461767845, 2.131737740833306, 1.2289466186044748, 1.059183785993618, 0.8954407393275675, 1.0866774642315318]),
array([3.8779001378782585]),
array([2.6578027031927043, 1.1299720246044815, 1.252608994830478, 1.4068758247148128, 1.047715849327966, 0.8617217822431481, 1.6880222692128628, 0.24302855485946429, 0.6044255249162527, 0.37953150861708773, 0.0]),
array([3.1694186902502737, 3.4215645146242704, 2.722347270486608, 2.4008562039152603, 2.2469243659245604, 1.447646186930721, 1.0674285694514558, 0.9124504437574393, 1.0539746648250545, 1.5126808917059285, 1.4731714575961774, 1.286627142742193, 1.7975528747958611, 1.3445785418915412, 1.0924739175648974, 1.1920591529548799, 0.5968638355392163, 0.41412928968576745, 0.6651078413945583, 0.5967837763199176, 0.4440202159058835, 0.31241139814217245, 0.6859740236349501, 0.46861015893197916, 0.0]),
array([1.277206520170704, 1.5192104555256973, 1.1980353505796117, 0.5162215913395174, 0.0]),
array([2.6810122040412723, 2.0006182465335365, 1.308303452822059, 1.6006455010268663, 0.18213774653256842, 0.6720557794289223, 0.32352632170568996, 0.0]),
array([3.0315370604569494, 2.768514285428521, 3.1292174280532965, 2.4060970137221447, 2.2639848604302126, 1.3119465990848218, 1.6837639997314597, 1.6698433738634642, 0.9969976737440222, 1.6992304115019177, 1.2907014458059916, 1.5502306079571258, 1.2668329092101518, 1.703865479664225, 1.0481825309259118, 1.3694464951820038, 1.71009863141028, 1.1453939072318748, 0.42307862209070435, 0.7531283411747272, 0.7333535888229431, 0.17848560138807157, 0.38058445620981013, 0.4393001894650403, 0.3769861471458077, 0.7589985292389986, 0.506777419719766, 0.09955689892467547, 0.0]),
array([1.3517575370504191, 0.8327612869592765, 1.4438242958310596, 0.1930889333682544, 0.0]),
array([1.2025154839439407, 1.3655017412581694, 1.1076756732272832, 1.4480009085975882, 1.667330992889962, 0.8432779411177159, 1.7512633239164812, 0.6385665430852566, 0.3714244614390895, 0.22485569084312662, 0.15439174409437717, 0.7221686301061909, 0.2723220106357579, 0.13105321232141398, 0.6910383271839192, 0.35164328747805035, 0.15970281624232896, 0.40959129542465444, 0.05509439726723328, 0.0]),
array([2.6524764926784843, 2.744533346224898, 2.5968993941717025, 1.6353275644714176, 1.2346322623172814, 1.3903342661606568, 1.714524364404206, 1.5391399168529436, 1.7738680579682258, 0.6294754716459656, 0.7122484859890699, 0.6711387025535175, 0.4157467818896757, 0.0]),
array([2.760097528674997, 2.2641689819352764, 2.4702090663124037, 1.4463881584104967, 0.9490307729242489, 0.4555439648274663, 0.5872520349880195, 0.7479490650201853, 0.6546856815591023]),
array([2.406504387052032, 2.0093886504218372, 1.1366326201892352, 0.9129884392434512, 0.22141356330091255, 0.0]),
array([1.7697996573911177]),
array([2.6884505121669133]),
array([2.5927931262589463, 2.068520244172157, 2.049514245784043, 1.3638615271661472, 1.0742674953085856, 1.7500260028267465, 1.375835639912222, 0.9602384209595799, 0.9063761132714596, 1.4684170147281164, 0.6472881068488568, 0.27823060409774814, 0.3405080631771516, 0.20509582535662818, 0.4628992210811044, 0.0]),
array([2.439207903129635, 2.3491332218501735, 1.1491051947939641, 1.5705383115983256, 0.8691914282953703, 1.5003716576856096, 0.9882276312030354, 0.9658352064451331, 1.5833246215225232, 0.9663887067733256, 1.0146948533982671, 0.9380500834455929, 0.9272381983749695, 0.7339430031956583, 0.765297733376187, 0.26135825927258804, 0.59013729334569, 0.5339388912554021, 0.24571012215078114, 0.21359510084532574, 0.6923585467995678, 0.2935375127291105, 0.3459275182146379, 0.31546524845276847, 0.0]),
array([2.505585135284462, 0.8607254718558902, 1.074518980310485, 0.9112843349452753, 1.3923803055922925, 1.007760118401442, 0.2196867040049144, 0.6299287978431936, 0.2922473254853208, 0.7167693045726811, 0.18968162693460522, 0.0]),
array([1.0692374792154302, 1.4652989475770875, 0.6284760004657055, 0.7319806928864648, 0.7262836905198017, 0.0]),
array([1.2687902391709085, 0.9840570278832571, 0.4312768394638924, 0.75694414037765, 0.0]),
array([1.3741569157723377, 1.6380234970813918, 1.5215075185858864, 0.7922666178112179, 0.9299731216049709, 1.6764165896100032, 1.6114959886770133, 0.5142437283409428, 0.29361866314977403, 0.35190654518665865, 0.6765638929575395, 0.16426944472499605, 0.3708581990849655, 0.2985495722512964, 0.0]),
array([1.2233309800462493, 0.0]),
array([1.412557564883515, 0.25622049870671393, 0.0]),
array([1.8243155128945574, 1.595411253681772, 0.824973252354628, 1.4634737903770652, 1.4327345889972214, 0.41952625249386144, 0.7703595153588734, 0.40787325487399434, 0.5855820094727764, 0.7783583444477674, 0.6219922051612565, 0.7292187852629235, 0.3608173596791596, 0.453022161459386, 0.705797005135062, 0.44639507071202816, 0.4349495834269212, 0.0]),
array([0.9092167611519224, 0.8487399085276958, 0.0]),
array([1.3453409583254943, 1.1386020525795701, 0.8340941870204347, 1.2430087014625097, 1.319760422899481, 1.3982656097221626, 1.4376665113022546, 0.2392899224774926, 0.547560409544004, 0.6845288958790147, 0.5233904915981034, 0.16112758943919203, 0.4696070500953573, 0.6474304902911824, 0.7260402010288218, 0.0]),
array([1.361822986047316, 0.0]),
array([1.0886304052486193, 0.3865292892194425, 0.21501801544772747, 0.0]),
array([0.8579037014473028, 1.3461647939438506, 1.5000893961151873, 1.0188526480316265, 0.17678157716952447, 0.7112668690429782, 0.6566470545861834, 0.0]),
array([1.2728370899605594, 0.419940500936875, 0.25514614709576544, 0.0]),
array([0.31837196321694916, 0.0]),
array([0.9011113628072238, 0.9167296847527814, 0.8147411700537005, 0.9546357245966061, 1.1968191028747404, 0.42520435786136374, 0.192240213474695, 0.34371437362394125, 0.7351265289213383, 0.6074900127326873, 0.454242801809942, 0.1943423358009273, 0.0]),
array([0.9897081049423643, 0.9756780273524142, 0.396875040028489, 0.6918160681477353, 0.5009997261310541, 0.6302955809510562, 0.0]),
array([0.8253780667095708, 0.728297908174314, 0.4081129852987369, 0.0]),
array([0.8349844289381012, 0.44434516944660724, 0.647903333162574, 0.5066492032553487, 0.0]),
array([0.5942059368286443, 0.7617992389582172, 0.40133154509565405, 0.3279113691547633, 0.1948888645637432, 0.5334477298394487, 0.2569771164606861, 0.0]),
array([0.2418747128554971, 0.0]),
array([0.4278836204116003, 0.5416588086395725, 0.3218733354323475, 0.08611575609712993, 0.0]),
array([0.5279603635851446, 0.0]),
array([0.027832286944965506, 0.03780882889489173, 0.0]),
array([0.07511787082459037, 0.0])
]
d = [data_1]
names = ["9"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T4', 'T5', 'T6', 'T8', 'T9', 'T11', 'T12', 'T13', 'T14', 'T16', 'T17', 'T18', 'T19', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T33', 'T34', 'T37', 'T38', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T48', 'T49', 'T50', 'T51', 'T52', 'T53', 'T54', 'T56', 'T58', 'T59', 'T61', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T71', 'T72', 'T74', 'T75', 'T78', 'T79', 'T80', 'T81', 'T84', 'T87', 'T89', 'T91', 'T92', 'T93', 'T97', 'T98', 'T100', 'T101', 'T104', 'T108', 'T109', 'T110', 'T111', 'T114', 'T115', 'T116', 'T118', 'T119', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T130', 'T131', 'T132', 'T133', 'T134', 'T135', 'T136', 'T137', 'T139', 'T140', 'T141', 'T142', 'T144', 'T145', 'T146', 'T148', 'T149', 'T151', 'T153', 'T154', 'T155', 'T156', 'T158', 'T159', 'T160', 'T161', 'T162', 'T164', 'T165', 'T166', 'T167', 'T168', 'T169', 'T170', 'T171', 'T173', 'T174', 'T176', 'T177', 'T181', 'T182', 'T185', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'T193', 'T194', 'T195', 'T196', 'T197', 'T198', 'T199', 'T200', 'T201', 'T202', 'T203', 'T204', 'T205', 'T206', 'T207', 'T208', 'T209', 'T210', 'T211', 'T213', 'T214', 'T215', 'T216', 'T217', 'T218', 'T221', 'T223']
def get_taxa_names(): return taxa_names