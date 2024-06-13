#!/usr/bin/env python
from numpy import *
data_1 = [
array([30.38481906784885, 28.606882907596052, 31.774190121597986, 30.97814058246311, 33.43687637012398, 33.32745212671527, 31.397425856414127, 33.41848583526938, 34.79125578053258, 32.457985079578705, 29.82537850545779, 33.08821573871094, 29.047342824432242, 30.152733545340894, 25.20250253345026, 27.150556258629177, 25.347077832195293, 23.945274535915296, 23.64253194470438, 27.510103793419958, 26.415556248771452, 25.983905112736327, 26.156073393857028, 25.97653718570725, 25.767863469538035, 24.264290001943447, 26.551034948524475, 26.35883865741112, 26.12145919337398, 27.33332238149105, 24.90400669395742, 26.070752663329607, 23.455106441151543, 23.935421512341534, 26.93298788727922, 26.871468160990567, 24.133888101108447, 25.261373136320344, 28.04823172211321, 24.906492923048013, 26.097791534393746, 24.118966209870024, 24.313485535169036, 23.423796996692566, 24.55462012568966, 23.471198564705457, 26.656800069388265, 26.38271145208903, 27.79985171425865, 23.70463392531632, 27.318116829376528, 25.304810073362454, 25.518379380643303, 23.973349149747023, 26.13289340529758, 27.376690112319633, 26.157779096789703, 27.856386527293573, 26.5715480471469, 24.9058862490535, 27.81028425957182, 21.97896439210514, 22.68376195781972, 22.83231055161567]),
array([27.100406811282657, 27.558449020177232]),
array([28.294141298877094]),
array([27.757918087355506]),
array([24.071199955781164, 27.343913648032135, 26.476769331775206, 23.585394962532266, 26.862314512921238, 24.54540794218189]),
array([27.515813278713892, 25.48386413830036, 24.918500022115868, 26.936246295911495, 26.405664453028407]),
array([26.968731294136276, 26.909122712644958, 27.567131490250755, 27.17973027479753, 26.290403148607368, 27.31312323587186]),
array([24.71481447003851, 26.740399686020375, 25.839656946816056, 25.641951647322355]),
array([24.93420556000106, 24.43915239244139, 24.12064476476326, 24.510325538589207, 24.074693881524663, 23.7541358831677]),
array([25.53683625953711, 23.67240938578105, 24.815681149555594, 25.78492833765278, 26.106165819627286, 23.906249695515555, 24.38172587790964, 22.524081501475727, 21.711839495209937]),
array([26.208098718335, 24.440908259057075, 26.29765056112202, 25.111904307176427]),
array([24.483763403387698, 24.21584012071857, 25.884879157809067, 26.268003793062537, 24.46727943919281, 22.76926993124145, 21.861109501614262, 22.419901357407873, 21.616693550658965, 22.671853655796983, 19.58117495873281, 20.051520201504236, 20.03408388860747, 20.417954527361495, 20.047550942822618]),
array([25.604298455018217, 25.613389231275416, 25.792539444470577]),
array([24.74340629926266, 24.464165177012546, 24.159818870297634, 24.911494245737497, 24.155264774027383]),
array([25.18099646422969, 25.397645289259817, 25.232928780895573]),
array([24.06532097321177, 23.71728547418111, 24.64455189663238]),
array([23.75371579456826, 25.079623264538892, 24.202318301216433, 24.375078079010365, 24.382781861267343, 24.962453703291303, 23.927781114360013, 21.02542736863672, 20.637066454564096, 22.20285962603207, 22.11190051348108, 19.64814550104777, 19.90284778793636, 16.793169397455568, 20.34249517612833, 19.863112309348406, 17.47029020247168, 17.904682949420987, 18.522240509161417, 16.830108290433675, 19.375575551999784, 19.51641086393183, 18.1438377928757, 19.96325991439832, 17.359820053377376, 16.260495931852937, 14.172730174572724, 15.363462369210936, 14.994597622689767, 15.363527241874527, 15.348011918438413, 14.734033825229151, 13.99409424617263, 13.871297163353853, 12.40208955585497, 12.971671826213624, 12.247084437997954, 12.781028720850298, 13.696435795829302, 7.841002182419406, 9.927897850161859, 8.876952627004334, 9.002476613176709, 10.742010217064719, 9.637483804953803, 11.619260522007997, 8.195998064814981, 9.633823408895985, 10.413068442118465, 9.90415635781743, 4.56094488310671, 3.0782751300871256, 2.585039178826934, 2.8900058256425623, 2.4110900087461204, 1.0792297877012667, 1.3402236812151342, 1.4632464688020812, 1.5874721152115623, 1.7522948562300387, 1.3515411833592563, 1.7971766068946287, 1.4402116558569653, 0.455774569779661, 0.0]),
array([23.81046790335794, 23.895549320379136, 23.642239960522353, 24.319428530754454, 23.380813212613162, 24.050514965040644, 23.311292190200387, 24.43383168741425, 23.743670086174003, 21.480959228471693, 21.11646050271436, 22.879071160313746, 21.282324474396468, 20.921097046502155, 16.812934348423273, 16.50099752599424, 20.11390937845691, 18.544566503204226, 16.403857611699635, 16.679104111691483, 17.78574971018171, 17.514750414791994, 19.250928478600386, 14.75356390051211, 15.918029912160438, 15.375789282380472, 14.5602228440893, 13.981659363627019, 15.65592560068292, 14.188616980848415, 14.079136628308417, 15.776955973875932, 14.024646833005606, 14.689954547999813, 14.150745299703395, 14.52159255478717, 14.510798333033263, 14.089129212544037, 13.734587137774815, 13.708839244673587, 13.076622957497507, 11.786303094954896, 13.099879739959917, 10.683540479368784, 11.259986473196468, 10.780076471125277, 10.060470588379424, 10.332132059750625, 9.328298692561827]),
array([24.8934123286796, 25.157284563989904, 23.188914451614288, 24.489656431855607, 25.243979439161343, 20.97737006835601, 21.894022420728707, 20.89902457258002, 21.61307138872584, 22.546382020140335, 22.28994158697736, 21.494813867715933, 20.124373457843998, 20.360613123538524, 19.366702790863116, 20.09610349404156, 19.258923126255283, 20.43187513607176]),
array([21.599090259098467, 20.95238417527257, 16.86791852948554, 18.089668769111796, 20.25216490209327, 17.47554102433678, 17.30001363165968, 14.329618077183168, 15.091135375399505, 8.347553873225035, 7.9982960725313355, 2.8955824322846198, 2.2978734660258087, 2.0588929706805428, 0.6493246623217206, 0.0]),
array([23.839302179652222, 23.452088692798434]),
array([23.28178019132884, 23.043900646796914, 23.100722549043333, 23.62346850790738, 23.86503837548231, 24.048890052772276, 21.214170554300953, 22.48882721985759, 22.020638720025907, 21.053691580684784, 21.8663799559429, 22.458423142234317, 20.649779291895328, 21.889347990041177, 22.245833604096642, 21.731403992592437, 22.265667369065095]),
array([23.5977434966845, 23.31692132892668, 22.32970606231055, 17.31271554454886, 17.3291891399557, 17.408599052360504, 17.00891974339456, 18.740928840175073, 15.744342901480762, 15.398466822758124]),
array([23.096831505526847, 16.159028041148574, 18.281909636028093, 19.242845023679, 12.980471937509167]),
array([23.65892403683338, 21.316632668047713, 18.186234505792726, 19.805487674882897, 19.439503426325086, 18.85613465004877]),
array([21.312296053507925, 22.14692344886083, 22.682369241963414]),
array([23.093880435653222, 21.797689945364297]),
array([21.75988702025266, 21.714861085112652, 21.973363685245257, 21.27413147275306]),
array([21.660591755873654, 16.87231737561099, 19.024577286835033, 20.28235424956114, 17.85610868313113, 17.793819709895708, 17.178845420286017, 19.273399559235305, 16.77136346675733, 16.98051656666882, 16.588418300536535, 17.060631690846506, 15.642698264015628, 15.502835838643206, 15.791504733643597, 15.549503069439554, 15.245602497171042]),
array([21.372512870806656, 20.730889777679355, 20.689773894417865, 19.369744568802822, 19.55895053972665, 19.076652562764437, 20.08072042405395, 19.38263453041169, 19.180547560921287, 19.20817306097393, 20.233909008016738, 19.940961741044845, 18.95692012012036, 19.208431459478447, 20.194050784553273, 19.671585906720026]),
array([21.140341392514095, 21.20046655335269, 20.82591122294627]),
array([19.102150088080226]),
array([20.500277174312217, 20.99951252709111, 16.597931379504644, 18.115661834644335, 19.986533761775377, 19.459502384406353, 20.33073499567244, 19.69917483463601, 20.166771925107756, 17.310974488767123, 17.30668150232214, 18.458379557104927, 16.502083770930113, 19.78609390462405, 19.572795747408765, 17.707013102421218, 15.553243648632497, 15.05440425952137, 14.661252955341705, 15.1415344099852, 14.346380250367565, 14.86194075778669, 15.271048263069583, 14.411725314439872, 14.778934930564343, 14.998987821284524, 15.623624220964121, 13.84817892357167, 15.59604764217973, 15.143080988654974, 14.71103339511722, 13.963440332580689, 12.408658887177022, 12.280673836892518, 12.38381617421988, 11.385822745788188]),
array([20.7343973144341, 19.61092246564184]),
array([20.58712435604567, 19.89382786577766, 16.386380652278486, 17.660369626792914, 17.7224542510787, 18.772190758135388, 17.730254213885477, 20.120357946070275, 15.32191915420427, 15.723610895860643, 14.450549958282128, 13.925446728726913, 14.998287853531115, 13.824283564982629, 13.545532622426414, 12.634036675702879, 13.251251585805585]),
array([19.093603344341876, 17.608802317315178, 17.763752899600497, 18.088674779353884, 20.202663674939235, 19.002320875222644, 16.03921334956223, 16.223937248754865, 16.77357124490807, 20.154984785390734, 20.111530888793055, 19.48937692293552, 16.739926700560552, 19.528458506720177, 19.432892254823518, 19.276020165034737, 17.162525222589615, 18.08171653642289, 16.593441910363175, 20.10009325438696, 20.209553120237977, 15.277819988526367, 14.076744171930388, 15.966889286793439, 14.508914683770874, 14.066835482606221, 15.368146302849507, 14.147350894366422, 14.681180669774193, 15.281583868620904, 15.877469629202219, 14.547490074151833, 15.162955298148054, 14.801065645365611, 15.542951875473118, 15.928391877494345, 11.772715219633389, 13.785747164262958, 13.567891880276653, 13.526040140311952, 12.226933913879884, 9.34864294226361, 10.558596071117663, 10.582775120405895, 10.18672721646608, 9.134197652197077, 9.812531091151714, 9.884291266329752, 10.028129145020296, 7.305411308763204, 9.778089866904951, 3.971634875459446, 5.1174421114349915, 4.257255211704475, 2.591916208903071, 2.958365315338765, 0.8130851517273494, 1.1351062722816256, 1.378381725586631, 0.8376317419773935, 0.7859236287854452, 1.0086075313466631, 0.0]),
array([19.62760187253001]),
array([19.82817505408008, 17.204314544016235, 17.682217607194577, 19.40261084307257, 19.17327756350821, 17.09689755997946, 15.761756336374619, 14.592025688087656, 15.696001981621064, 15.330286634514076, 14.69734546084951, 15.382986080459004, 13.585592060372031, 12.874727393242955, 12.695329982583525, 7.603273796347956, 7.2768819594675564, 8.906893564978425, 7.347368941815669, 7.246966902025633, 4.671765295991322, 3.597099170106402, 1.489979711011314, 0.8973073492504043, 0.9940620999464039, 1.138350026004139, 0.7714873803438527, 0.5978550831000857, 0.0]),
array([20.168782668277775, 18.394233259226674, 20.124630178465498, 19.90283827607497, 18.994091961465116]),
array([19.47705554235401, 15.979205430897922, 18.215326921382115, 19.656946650552293, 19.276517365973206, 17.7826531160597, 16.48363133080837, 18.91209400064055, 18.34845794060364, 17.67688753409577, 17.353265576615517, 17.065836773903886, 18.643973505789965, 18.462273590708797, 16.47053281361488, 16.341907257107035, 18.277097305759394, 15.983662916115382, 19.236304049395788, 16.66976583380256, 16.330662507565485, 17.69393356438206, 17.18998874957002, 14.953596874160029, 15.083129243199132, 14.298836019719847, 15.102673127653485, 15.361431449574756, 13.850161686140673, 15.598461732713812, 14.59328263226385, 14.629027499381923, 15.218382592410348, 14.37113589838252, 14.349416006744509, 14.945821844472395, 15.275414325470276, 15.726509643341995, 14.09970667770501, 13.85473166462603, 15.367015152986488, 15.281773280477635, 14.581559595167867, 14.295523219059037, 13.156988949895794, 13.43980141326965, 12.682427235228348, 13.236304648169076, 13.593682938418372, 13.100553509461605, 13.674681503013849]),
array([19.673214746924796, 16.590013921437468, 15.598048729433888]),
array([17.215168668920914, 17.6268001518945, 17.61575890771193, 14.319067003984001, 14.507107266480867, 14.632704900031229, 12.353950154590814]),
array([19.179049562471388, 16.553309730665397, 17.76874329875162, 17.46776021069325, 19.14028701603544, 17.390550396409672, 18.879485855787067, 16.877736066764434, 15.938981030601411, 14.917044599045706, 14.446633290971734, 15.319136517557011, 15.43915262412132, 15.082486467748835, 14.799791608380268, 15.946298261525655, 14.354125351141953, 14.446398168794095, 15.48670991102472, 15.43552926595742, 15.594854498736156, 13.440354337269184, 11.8700037976698, 8.057034591591599, 8.495479677035096, 10.660691456721953, 10.824826737260784, 11.397429145162839, 10.67848451617729, 6.754930260083284, 3.846465884107311, 2.9102134060276548, 2.958618632381991, 2.4534016731349713, 1.3222916578357204, 1.5254156550726121, 1.4965856198652583, 0.24973648656243752, 0.0]),
array([19.671285947417285, 19.395398114012334, 19.148265396005804, 18.60729805191176, 17.26855380317026, 19.73061222177809, 18.165575319842322, 19.647383608154712, 18.668087004455987]),
array([18.566563772704452, 18.972652407711337, 18.017617407618477, 18.434282410872896, 18.67654880460359, 17.616592691077752]),
array([19.001034124021512]),
array([17.66114430996438, 18.723251372857977, 17.636830847873142, 18.852782951081277, 17.69364157185971, 18.99777980549306, 19.079085437797957]),
array([18.596320798741438, 15.084568142831683, 13.906312738659418]),
array([18.540277671458483, 18.654163972199903]),
array([17.889584440179, 16.835695520154236, 18.54182581750658, 17.188302763186872, 17.559101873252917, 18.19240057060335, 16.20390046457843, 17.665001338563588, 17.847950296435826, 16.779564004434643, 16.352993100001004, 18.015712992107083, 15.130154271295046, 13.979161633485704, 15.580311079942923, 15.455342555007068, 14.880091128685068, 15.645190189671931, 15.726369000863263, 14.52948175663963, 15.23385571544153, 14.46764670985591, 15.48141931540348, 13.648293186238966, 13.141452776575546, 13.347365747753482, 12.572429691241798, 12.932441400927903, 12.8294770327038, 10.985628563135169, 11.096607810883587]),
array([18.26555358392124, 16.07511755246341, 15.334402380797709]),
array([17.45865503105505, 17.179023334526853, 17.381955965058058, 16.41732390185717]),
array([17.437649504284938, 16.594350401771703, 17.19651541717781, 18.335200580630335, 14.050395792753452, 14.723747124901777, 15.706374137016498, 15.065578031570846, 15.754840262012245, 12.839205550870382, 12.479884969530545, 8.320991649490724, 10.310384302578433, 8.909494833568866, 10.047710423546075, 10.660631101444697, 2.9595742106262506, 2.7422518293706055]),
array([17.181953696794892]),
array([16.31223094429858, 17.08125280077556, 17.409085306073777, 18.279229485245846, 16.67793121542443, 17.40953757898129, 17.573676962055934, 17.463644321645347, 17.350045896818422]),
array([17.5243171621574, 15.827590112355653, 14.757053376757902, 15.72975621091422, 15.509108731534726, 15.813915326489418, 15.730755279275336, 13.24628913689687, 13.11491601521501, 6.766109366298803, 4.644428905508211, 5.0032752699741065, 2.848190027031994, 1.2721284831304824, 1.2857275895236402, 1.3389077710950585, 0.004364099063685012, 0.04959306436136901, 0.0]),
array([17.88998593665997]),
array([17.01663169335654, 17.819545332764292, 17.47202438877957, 17.060025840906405, 16.193539517185187, 17.079425833218995, 16.831567946291237, 16.17643005448033, 16.346701700342596, 16.528145201917063, 16.708253233079695, 15.005883095912663, 15.613832003334107, 15.294046129802416, 15.03440025644284, 15.55265369901208, 15.679996909004833, 15.652222360006427, 13.94626658437634, 15.220277255035882, 14.990435834111357, 14.351170415250973, 14.277575226132017, 14.474194156635932, 14.052316446795231, 14.82478596192125, 15.613111420868977, 14.65371717107434, 14.997315900772794, 14.186356538713751, 15.447634497924431, 14.514210538594597, 15.948342726713564, 15.146415373015236, 15.225884723078481, 15.767118780943704, 14.20512532624843, 13.254250024860236, 13.502655599932531, 12.452918042204345, 11.774502099946046, 11.830575725791459, 11.706619069492369, 13.578914530811959, 12.250801516913523, 13.656590903947393]),
array([17.262406263846746, 17.69613556353023]),
array([17.24134655626034, 17.358151878641973, 17.302088542717005, 16.53597709004059, 16.545534619657257, 17.385190352200993, 16.207651336876157, 14.192232346844735, 15.431669640746389, 14.440123460277793, 14.030515498051894, 14.722815700917145, 15.177155818446408, 14.503051907070956, 15.11972166624896, 14.741383589313468, 14.548011310679032, 15.340288165339452, 14.546048043699495, 14.8781032870151, 14.694815418040946, 15.348646926687854, 14.36649483594108, 15.64897913742825, 14.9362953955793, 14.826805118170954, 15.952388933621494, 13.42051012524433, 13.40639990306778]),
array([16.755105817959876, 16.622338391572598, 17.2592732942788, 16.183877615295156, 17.250308274063933, 17.07092781888769]),
array([16.488791740652893, 14.255874735139958, 14.777868105299754, 15.800410542425317, 14.14493042417095, 14.443030694809002, 14.00333884805304, 14.503849771967527, 14.887501068895155, 13.062727538949693, 12.260493808370821, 13.08768056745795, 12.964945754609175, 13.285669709464484]),
array([16.45225683295415, 16.418105759037065, 14.270632621857274, 15.199282788331686, 15.033459370456141, 15.204361390660075, 15.268976006917882, 13.941741164357824, 14.881174916970837, 8.228304523257236, 10.346934960439842, 9.355201341820278, 10.529044364944912, 8.787461673083683, 6.161849473974584, 5.961514828958536, 4.734061688126584, 2.5856744040110935, 2.7962794241457636, 2.600420456010184, 2.1678703448367544, 0.9893637356035849, 0.0]),
array([15.624902093153409]),
array([15.156141585586331, 15.265470883772274, 15.47347611138912, 12.02729003074537, 12.585305198658492, 12.723685473228397, 10.866218450136373]),
array([16.013491326841166, 16.62621832199818, 16.190168907163276, 16.357747799187734, 13.92934744991526, 14.395008121758545, 14.935811107489064, 14.685326499151051, 14.16694185032854, 15.174609722801913, 15.363455329642711, 14.94393332446732, 14.980561008946719, 14.183585752457123, 14.046875922377176, 14.123139953596619, 14.477241541823263, 14.83497164899915, 15.59447148908999, 13.967699134080323, 14.133493253775345, 14.838853311192787, 14.348553638874698, 14.036869144690955, 15.315293717096162, 12.889666529217845, 12.170218833631427, 13.661967861055585, 12.43504586397376, 13.660138906716679, 13.344473727202754, 11.49009365335788, 10.232131056818389, 10.527533269095875, 10.178312508848693, 11.455908560630984]),
array([16.078901493423885, 15.396219471290863]),
array([16.129594320482692, 16.059948437619664, 14.934229144050901, 15.728720692691859]),
array([16.03550648741447, 15.891935733491174, 15.476383202549854]),
array([14.08679449446917, 15.876051303318365, 14.722934710890064, 14.582086053027039, 14.42874811378206, 13.74748440143271, 8.221634104092068, 11.418207808841096, 7.4373288029088, 7.5122734818348125, 4.4218000101905535, 5.191328321313779, 3.156049058451429, 3.348855899474683, 2.292425824648763, 1.9873958748789367, 1.5458453897698803, 1.6743729483224152, 0.0]),
array([15.989763839845248, 14.645027668905117, 14.760617094899564, 14.66387019642474, 14.637728803017444, 15.77156662862443, 15.595965303888171, 15.354814520346332, 13.97498649656463, 14.965911814364837, 14.650637112412477, 13.887897882586547, 15.754067875673519, 15.556121211130858, 15.36502894041588, 13.993936944609345, 14.426255850833378, 15.689308984851017, 15.817427053778138, 15.604591026191853, 14.001470046189835, 15.064704501142252, 14.100061870614667, 14.149772565785701, 13.892094308531197, 14.804128183015084, 15.382232112031653, 14.863991618446871, 14.049503717690156, 13.93184651378344, 15.114940787882562, 14.922669998491845, 14.254570307674207, 15.116513570077016, 15.454708533222602, 15.369352440453842, 14.9377446300045]),
array([15.86512264181112]),
array([14.276623676413083, 14.465609702201355, 13.950634998769674, 14.311200619513086, 14.880577042352384, 15.405668985658197, 15.303504788081295, 14.039483816136098, 14.526949390971447, 15.552646350189319, 14.79032777689446, 14.273509417847464, 14.19872037255032, 14.13156025061867, 13.586234036908586, 13.395403579223967]),
array([13.975713238605696, 8.118882908132967, 10.258552291741802, 3.8356936009940727, 3.037702922471026, 1.6253303879329677]),
array([15.468007523353874, 15.436566855540255, 15.53753077361741, 15.437528802579182, 15.462513215144579, 15.513767096670334, 15.50700693667324]),
array([14.322680705867924, 15.151453664917797, 14.95057978557958, 14.891638875069576, 14.231179898715835, 15.208641551301643, 15.403752517396777, 15.22455610132659, 13.690056828462048]),
array([14.311580227207866, 15.443687833305091, 14.582073937337404, 14.667861275554001, 14.890158245065788, 14.93608072911407, 14.608872181859963, 15.42596498921285, 14.977900123560698, 14.817851871241247, 15.287910617525563, 14.632893832322635, 15.000887164744489, 15.41783228966842, 15.31096335097512, 14.969377744997118, 14.398269991157816, 14.16158462687108, 14.698564223545485, 15.450121499558659, 14.492052525029276, 13.98190174118146, 15.441889601212413, 14.935427727886125, 15.501174382717531, 14.542822384465678, 15.351669713005244, 13.920998292569497, 14.833371975782542, 14.467479132002513, 14.967923007015456, 14.955566234990131, 15.510267923364355, 14.048029153500053, 15.068772051262938, 15.105305658826795, 14.914506206140144, 15.215157884020767, 14.995922561992368, 13.572730839947074, 12.666850372864353, 13.774305567245472, 13.052051391169764, 12.84476091061445, 13.65021780521461, 11.877542884698183, 12.979156220385422, 12.534060959839799, 13.573358482302071, 11.699334917697211, 13.421898516193156, 11.513781726906087]),
array([15.16884235262006, 14.697029488006859, 14.764132619208407]),
array([14.910751359636429, 14.06443559824474, 7.281784949710568, 4.203395797471596, 0.0]),
array([15.176235966645182, 4.7889542983515625, 0.0]),
array([14.242768237939453, 13.953759787043696, 14.587443012259985, 14.642270818115371, 14.02176075766751, 12.79895488184204, 8.486451732615702, 8.734259218203, 11.497687161211566, 4.277995001414749, 1.1682602531301196, 1.7674226580541368, 0.0]),
array([13.626711892492045, 13.475839464259733]),
array([14.080024934732085, 13.839390294670062, 14.450911198392657, 14.458674429910124, 14.52846080354642, 13.522663029890662, 12.975556127612606, 13.627666817321154, 12.797666782156917, 12.582161988184044, 13.402229396359091, 12.138828333119926, 13.653367541365474, 12.561034987690883, 11.582395039248192, 10.806619732389393, 11.020556516928837, 10.987559380155774, 11.378896456440756, 11.230354083322512, 10.464290790921092]),
array([14.151099635537232, 14.28661362987312, 14.009160177080139, 14.267925193321398, 14.087422705350745, 14.057616257238923]),
array([12.046093432361198, 10.192075473852153, 7.804780813719006, 10.972810159960444, 9.53478289221219, 9.115185050576411]),
array([13.268720849937322]),
array([13.730612599614286, 13.581866849321525, 13.642960491360892, 13.46753405344592, 12.6296302788942, 11.71594807080521, 12.961851268164066, 9.036566661266395, 8.585429364361923, 8.771353525867978, 8.538182045694825]),
array([13.321892959552985]),
array([13.053028057788985, 11.683323875039815, 12.833661378083484, 13.199493248361405, 13.206896825519715, 11.401877981578394, 11.317867805866062, 11.116220733038853, 10.865696622475433, 8.857836019937917, 9.553791768566605, 8.901245232675164, 9.317196435832265, 9.469594710553668]),
array([13.171218979818923, 13.228325789840854, 10.815928248025964, 7.858793624198351, 8.241897014364568, 8.315662733903217, 7.545147617877085, 9.89633588581814, 3.6281689634753307, 2.693869266915163, 2.0776700054338333, 1.4823623758150668, 1.626203231283066, 1.7725811742494586, 0.0]),
array([11.742239456596693, 9.78043544239908, 9.56472546278114, 10.09685172126402, 9.388842539528966, 8.884637748579287, 9.05475828575029, 10.856525811490377, 8.450698584945798, 11.320498430341772, 10.721566180458371, 6.768222009960946, 4.287145576099764, 4.627169595326414, 4.151099716132922, 2.7231672231146553, 2.9166550966952887, 2.463583697546133, 0.8387100554751921, 1.4291208960914537, 0.8818876886815117, 1.3093091453603205, 0.8114988869815122, 0.4180466998317697, 0.27656101294003177, 0.0]),
array([12.003628901862474, 12.783500325761068, 12.16717201522616, 12.353629573421768, 12.882652082466237, 11.646524553443651, 12.660512844335132, 12.206016825231963, 12.493480499734996, 12.321754686338904, 10.352683048067334, 8.967734606836986, 10.818401116465354, 8.965562598293578, 9.760497826235053, 9.378842928326717, 10.675777220572225, 8.50318715383717, 9.533907508574465, 11.3708986503717, 9.033558453459639, 10.125348034547505, 8.941994591565017]),
array([11.642120385854978, 11.939776567264458, 8.558637503736222, 7.39408536665212, 7.644154877548836, 6.19238343348763, 7.1889323878177, 5.382865160086029, 4.02166678682479, 2.3390107763142476, 0.0]),
array([8.717489996854702]),
array([12.46017059017253, 11.952762192087219, 8.274366027437853, 11.366728133431302, 8.145208500144754, 8.4066141083126, 7.07318157518993, 6.3036726130682945]),
array([12.144041588876451, 8.017173667091527, 8.738044250402439]),
array([11.980367462649598, 10.470150602237554, 10.166574736592842, 10.497515314538502]),
array([11.002495469847624, 10.356251287785339, 10.987492756636158, 11.498694183377866, 11.570690576871721, 11.150432777023099, 11.036813562892679]),
array([9.191218824816787, 11.324700030023607, 9.217076709444337, 7.427794239512483, 8.334559754550476, 8.546296907760281, 9.371268182182439, 7.752994684163923, 11.23148889056562, 8.221624851370182, 8.503024206055713, 9.971495920411925, 9.479495401336765, 10.963665359565463, 11.602646970391751, 10.15803645327258, 5.335674230969388, 6.629346250981165, 6.6958676681749285, 5.831329191118814, 6.353674828665417, 4.1420884690430215, 4.964837716668827, 3.1883408241504054, 3.3193229969096527, 3.574364225139044, 3.174165056952283, 1.94981907884629, 1.888034052942935, 2.181355997015594, 1.9348403419395117, 2.415726738797214, 1.658262076700642, 1.4218006982876272, 1.546683915970792, 1.0368021383249204, 1.5523127854139098, 1.2297205770146675, 1.6073832578304825, 1.7319534717772547, 0.7293473924317935]),
array([10.974108817278001, 10.620194344146315, 10.959014253729418, 9.9046821240248]),
array([11.179183457353206, 10.79860764338296, 10.955184182181657, 11.136108758783276]),
array([7.721252864106914, 8.959525080871206, 10.709551265774044, 7.757291081575351, 10.035384075862384, 6.735251331629141]),
array([8.992315201860567]),
array([8.438553576500617, 1.0883131149016831, 1.1827323709540523]),
array([8.10464082818918, 8.196592133445497, 7.810878274796224, 8.01570227613744, 9.967166810676348, 8.731468194403302, 8.94404575130466, 6.148394553818769, 5.293794508964859, 3.0867060941470132, 2.5982717386059977, 2.3120050791030007]),
array([10.023089161992743, 5.661310988000897]),
array([9.160705431591389, 8.590676799280015, 9.142053709741345, 9.22909226102216, 9.089319104493988, 8.942192500817326]),
array([8.76147707232527, 8.636891948981965]),
array([9.090225330913531, 9.156063440341386, 9.020255931931967, 7.595910570290547, 7.503943513578928, 8.953293551082124, 7.678061948609856, 8.428001579779487, 7.355194254973279, 8.361364096159912, 7.494643467941226, 8.40165226541558, 5.9638872948802]),
array([8.758602263614527, 7.596484996929824, 8.786006208215051, 7.737307087276686, 8.336722926643759, 8.81797633996224, 6.785986730803538, 4.781787201649892, 4.577958780589638, 3.8918071869151043, 5.177162487427035, 3.8705122598476174, 5.227095441330461]),
array([8.290418601465543, 7.433241118049392, 8.647020606186022, 9.067869878490633, 7.437684413413993, 5.9664741812036715, 5.464347186151566, 6.868285452948861, 4.786923633921843, 4.631018763694978, 3.0211410073767024, 3.0098579797334453, 3.124985138037266, 3.585167212882108, 2.0781032976613716, 2.5583575820026407, 2.1669399306363757, 1.9148569917089224, 2.4161567988544546, 0.9910681561358692, 0.7839331021839324, 0.9330903847487553, 1.728852553237425, 1.5119001560625767]),
array([8.983857057758062]),
array([8.494969694741057, 8.090980471385887, 7.603511105072271, 7.282473855452377, 8.56658844468077, 8.492459422659191, 6.329606069213256]),
array([7.251841550246221, 7.9079394837857935, 6.2888050904845745, 3.718159164910123, 3.6881376169357, 5.259615071893261, 3.30759405020275, 2.5708113619271526, 2.074818358199928, 0.2870412453762377, 0.0]),
array([7.4456804720958205, 6.824110486304406]),
array([7.933773931923491]),
array([7.909141897805714]),
array([7.345236960344532, 6.982123973440508, 6.799864171910582, 6.579959088637937, 5.637321189740364, 6.509622670312887, 4.991012497948397, 4.554583909438283, 2.929414302927321, 2.8776583021970437, 2.064731139829222, 1.3787861219791997, 0.8330941695554316, 1.7788584069268925, 1.0473358753537458, 1.732995072586673, 1.6242497703365724, 1.5495280994081264, 1.346809494546083, 0.8785137811000792, 0.13489563086376444, 0.7638487536346373, 0.12382617662905047, 0.0]),
array([7.860031927664215, 7.81459386825534, 7.73483499678063, 6.5994910715824355, 5.825059146624611, 6.980200298269548, 4.8587340597677064, 3.680682148570486, 4.2709451687607345, 4.76505844519603, 4.292496161620255, 2.630521325261203, 3.1087816432661013, 2.39072381126868, 2.3253935399961687, 2.31167337072961, 2.5150020643505635, 0.9255068202379885, 1.477638842610255, 0.7966823890344645, 1.7376499525889315, 1.7399170397806611, 1.3098447896021779, 1.453235599925804, 0.9137392530232575, 0.6935705721602212, 0.0]),
array([7.66757890549004, 6.632901637784204, 6.190634087353674, 4.163217758571278, 3.8480480865191984, 3.7945760998595777, 2.8979433525219704, 2.7874205728710675, 2.8264608375381757, 3.4717367620587583, 1.99956829307976, 2.4127994529162535, 1.977931596410249, 2.07192137936031, 1.7741035645225949, 0.899734013163782, 1.5120749487257612, 0.9948053329174575, 0.22144154103989278, 0.0]),
array([4.622325694419638]),
array([6.582357705030311, 2.6917132347128714, 2.499101926173932, 1.7833590122758936, 0.8677563379227511, 0.0]),
array([6.69763619524029, 6.694807472360952, 7.058374212732068]),
array([6.742518060843601, 5.681631567000043, 4.837060620007151, 4.697226103959577]),
array([0.5311123322096014, 0.5703735122346495, 0.0]),
array([5.697060015588033, 4.406560965173213, 4.064303292006997, 4.484269175582177, 2.8318802544168835, 2.9171734051002574, 3.453080082379055, 2.6378750178293178, 2.906469665809581, 3.405757322458424, 1.9965892084768089, 0.9574070528017935, 1.205690410361093, 0.22451342422782428, 0.20418680599964112, 0.03590873582452557, 0.0]),
array([2.4678488919227606, 0.0]),
array([2.0665507197132156, 2.5194694771209725, 0.9519573174838273, 0.0]),
array([4.668019135657669, 4.978150395567427, 3.002903805625723, 3.448825240353292, 3.5439096322646377, 2.874094274289782, 2.127077394317121, 2.233502454396706, 0.9829822254046124, 0.9288037908438935, 1.1033957134444135, 0.811583522632787, 1.2115646156910214, 1.0890656025050607, 1.315463389841522, 0.13082062113892123, 0.6168708680974225, 0.40179330679937825, 0.0]),
array([2.1711635922649273, 1.9960793460301876]),
array([4.981386102682995, 4.626763196012471, 4.58083239531736, 2.847963065630016, 3.0743292446965746, 2.0041329057498345, 0.0]),
array([0.7933492870575154, 0.0]),
array([4.0128166025178205, 3.916300507917658, 3.8776502740832384, 3.8216128768344957, 3.447271307595307, 3.26654583888695, 3.493155995027744, 3.36872103264908, 3.5328022338747957, 2.923551238042015, 2.8060753193672188, 2.544513669516607, 2.452633113452128, 1.8140300357721824, 2.237329207077743, 1.939099876247199, 2.1504493296060407, 1.7772131057026443, 1.0485626750542285, 1.5529822625078078, 1.1139610712608863, 1.6203260781072675, 1.4511379380338716, 1.1001528867455836, 1.3833755489327317, 1.122969274674027, 1.1416831985110165, 0.7162123712490499, 0.45074109582086075, 0.23192191999282274, 0.5630023659992558, 0.4453690052936236, 0.0]),
array([3.0723236731567294, 2.926190670251271, 3.4466738897668643, 1.976862872555055, 1.6835770781965487, 1.0930524902030938, 0.0]),
array([2.524780235515715, 0.5271635609951607]),
array([3.2367857104945434, 3.362331854102989]),
array([3.4729749483276366, 1.516034930502776, 0.0]),
array([2.188767862368682, 1.9931545328386735, 1.642840118924107, 1.6100875287106837]),
array([3.036580872923354, 3.110480041158303, 2.7093044973154217, 2.2560932566888794, 2.022087899752405, 1.9317975607503115, 1.4133412464558157, 1.5310972140763588, 1.6346583558375014, 1.1613833966059701, 1.030439873457972, 0.0]),
array([2.7777358989083787, 1.944816792336574, 2.4197809313545204, 1.8723124322404532, 1.6051746539731504, 1.185090899560298, 1.6661732540091903, 1.1436142509980103, 0.7341820276703237, 0.0]),
array([3.2571176234531944, 2.843965922798925, 3.289984417044821, 2.12584960950606, 2.436963536936682, 1.9617386507147807, 1.4958843802261403, 1.434713056829353, 1.1162402605004045, 1.3272486069182465, 1.3419677032809036, 1.4907920997724142, 0.971385660933342, 0.5922105052926696, 0.0]),
array([2.8875560837634353]),
array([2.655513057062791, 2.9533745458464016, 2.7422545227494903, 2.853111903652691]),
array([1.7453310875068289]),
array([2.3816007822472045, 1.1392312337902106, 0.6006814684393949, 0.33127788922853535]),
array([2.8720294280299243, 2.617528105929884, 2.3203919734096714, 2.255292637901779, 2.1857747240640535, 2.4455171173788166, 1.9169917725946002, 1.9405275696150248, 1.7571776673578992, 1.1444868486847875, 1.1042997224866797, 1.6680649063927422, 1.48227962030475, 1.7224234457643246]),
array([0.32534317028542875, 0.37324117640757554, 0.0]),
array([2.6080282607870013, 2.3421498282809066, 1.4275333731007087, 1.0186679545436979, 1.2657458518547213, 1.6762620128592973]),
array([2.6166465003286508, 1.8231864177276276, 2.534022941282125, 2.48468125657296, 1.9824989015145653, 2.39081303611498, 2.4172914070802514, 2.5625317163376686, 2.323369566942546, 2.521472115271984, 2.082562391557618, 2.57512944732582, 2.23945039833829, 0.8120744136972634, 1.2320020603826252, 1.7242360253699585, 1.4274326271617186, 1.6514621657108366, 0.9000637998644786, 1.4010028867018904, 1.290822597922575, 0.9437157183222805, 1.2215267682967466, 1.7505635784963285, 1.724151670359253, 1.6231025221630717, 1.112569561512686, 1.170915151928223, 1.435504764374425, 1.480955378167199, 1.7183062996639176, 1.7290136327611176, 1.1528619288675985, 0.623641113894357, 0.22177579659785474, 0.0]),
array([1.8686716994512071, 0.0]),
array([0.9897585407833907, 0.0]),
array([1.965436514621963, 2.0547422543008715, 2.138816860753175, 1.4534164346866572, 1.324751948961954, 1.4912171223029302]),
array([1.81736035620492, 2.1730267883331744, 1.99009975032046, 2.039944532488955, 1.887847401295191, 1.7446022064762996, 1.7976386628667347, 0.9834111181175931, 0.9163981293028072, 0.7975971138049458, 1.486668719731648, 1.0523694376090018, 1.0457566431279643, 0.632943996260595, 0.27668026072215635, 0.6021817952133117, 0.10313004159378822, 0.0]),
array([1.943388961010044, 1.8215326194528823, 1.8168320776013671, 1.895465326214705, 1.9832266921654358, 1.1217764144826998, 0.9637250351043747, 0.8263342600383012, 1.5944740762888552, 1.0578862707075696, 1.37879821325697, 0.9570053386583336, 0.8702093341408901, 0.8407433127924098, 1.0288478831490793, 0.8629807328451102, 1.6797182427606867, 0.16026550882192325, 0.7459089223370463, 0.7076296836778221, 0.0]),
array([1.728631689961469, 1.483197571784479, 1.2526511000416567, 1.186709777394499, 1.1638681648440337, 1.7000146479012865, 1.2801367091446734, 1.3788597354870489, 1.6729836822298725, 1.4394329931621794]),
array([1.0101382705544237, 1.6316596402890289, 0.8880037224679435, 1.6325979343711572, 1.716706850998051, 0.19601007201942544, 0.6726670340363005, 0.3969494918238303, 0.5458025287691427, 0.0]),
array([1.8419498774412788, 1.4297832792597052, 0.0]),
array([1.8048680922279052, 1.0068971149765367, 1.0101607162474782, 0.9186146687071638, 1.6169938032374358, 0.40258605535402453, 0.379044028932603, 0.07618901629639699, 0.0]),
array([1.460587858889497, 0.0]),
array([0.470420824950112, 0.0]),
array([1.407539559301037, 1.4064790824077573, 0.45546979514659625, 0.0]),
array([1.219273021835213, 0.7847793896015455]),
array([1.0226679999103512, 1.2753236634680172, 0.8109575588740935, 1.1088108130978833, 0.0]),
array([1.4336686236773155, 1.4101538926288806, 1.083820525907629]),
array([1.3650387534085557]),
array([1.1676335934951536, 0.2979981319870966, 0.0]),
array([0.6613164076342181, 0.0]),
array([1.0375689355529618, 0.9460246210316978, 0.8091468973946777, 1.1539402211663814, 0.32220453683798966, 0.7605292917383245, 0.26451949812987785, 0.0]),
array([0.5151203805086753, 0.0]),
array([0.8047801282996812, 0.33377161433343694, 0.0]),
array([0.928574154753721, 0.9594925085209306, 0.9652438849595222, 0.0]),
array([0.9283090488360285, 0.0]),
array([0.6833869554041743, 0.0]),
array([0.1504054271330535]),
array([0.1605188740239183, 0.0]),
array([0.171266255640323, 0.0]),
array([0.17827738067688284, 0.1968275725861629, 0.05683640670413609, 0.0])
]
d = [data_1]
names = ["7"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T4', 'T5', 'T6', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T43', 'T44', 'T45', 'T46', 'T47', 'T49', 'T50', 'T51', 'T52', 'T53', 'T54', 'T56', 'T57', 'T58', 'T59', 'T60', 'T62', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T71', 'T73', 'T74', 'T79', 'T80', 'T81', 'T83', 'T85', 'T87', 'T88', 'T89', 'T90', 'T92', 'T93', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T100', 'T101', 'T103', 'T104', 'T105', 'T108', 'T109', 'T111', 'T112', 'T113', 'T114', 'T115', 'T118', 'T119', 'T120', 'T121', 'T122', 'T123', 'T125', 'T127', 'T129', 'T130', 'T131', 'T134', 'T135', 'T136', 'T140', 'T141', 'T143', 'T144', 'T145', 'T146', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'T158', 'T159', 'T160', 'T161', 'T162', 'T164', 'T166', 'T168', 'T169', 'T171', 'T172', 'T174', 'T177', 'T178', 'T179', 'T180', 'T182', 'T183', 'T184', 'T185', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T193', 'T194', 'T195', 'T196', 'T198', 'T199', 'T200', 'T201', 'T202', 'T203', 'T205', 'T206', 'T207', 'T208', 'T209', 'T210', 'T211', 'T212', 'T213', 'T214', 'T216', 'T217', 'T218', 'T219', 'T220', 'T221', 'T223', 'T226', 'T227', 'T230', 'T232', 'T234']
def get_taxa_names(): return taxa_names