#!/usr/bin/env python
from numpy import *
data_1 = [
array([34.42056612904309, 32.61313687673921, 33.99798151332537]),
array([29.259038309732645, 30.30345059282601]),
array([32.39873335444033, 32.42237514875991]),
array([32.534190691256704, 28.229117946233842, 29.007150534307016, 28.951472702605017, 26.67174443474932, 24.089084548731694, 20.974343552721404, 20.49386193785229, 21.559161498944523, 22.447354752041505, 14.472253611018619, 14.387018468913915, 9.79420080128299, 10.332378376760419, 5.4111456711873105, 4.315804223015119, 4.350310398300666, 2.3578128803427783]),
array([30.633423242693755, 32.56622733104399]),
array([29.01019496929536, 29.89500137828922, 28.848322918587527, 28.36378384735808, 27.99556164049086, 27.61705985403252, 27.881326013365893]),
array([24.8361447948246, 26.06589089325733, 26.461044757106663, 27.94587949265565, 24.32475118037025, 24.792850335715304, 26.989016436924732, 25.9195072592069, 22.81647997159976, 20.98513594053625, 21.177231040890852, 22.221788412871497, 20.926340730851273, 22.192021215427303, 21.876051256836508, 22.347498458925706, 21.536093536598234, 19.266048347799156, 15.218171840224674, 14.214024867171457, 15.392130818630234, 14.669922075433178, 15.319711702569988]),
array([25.16387521608032, 27.27666619708597, 25.950363372722435, 28.03190465345617, 23.628892221105453, 25.31303455054771, 22.299536114830378, 21.949668418649644, 21.352169988077463, 20.821094483648675, 20.675271372194565, 21.00029522273484, 16.238743062043866, 14.7910258788101, 15.19546205359028, 13.998093596647166, 14.549894344653024, 14.39184297717693, 14.559382610803857, 5.627262692958073, 5.437124260029075, 5.922118944024019, 6.076158706890121, 6.616555834226534, 5.2325464444075305, 3.6920339045171566, 2.5536294343666155, 2.039918643790744, 2.526741474272383, 1.667178479434444, 0.0]),
array([26.04556078385382, 26.420768314383977]),
array([26.295620335646245, 26.14248421396444, 23.748117638754884, 25.684029077901908, 24.445219299527572, 22.19929451249027, 21.033305329458337, 22.25556265918229, 20.815859315195766, 21.686855294953173, 22.35279129445002, 19.80720900031333]),
array([25.009291545596852]),
array([25.3516350720818, 25.431414724909295, 25.395433159120877, 25.504571591255335, 25.654836631281356]),
array([24.684682912407673]),
array([25.074777514442857, 24.89892086331999]),
array([22.106011616238735]),
array([24.249891070082096, 24.15875125324931, 23.1859189914186, 23.397410726572204, 24.332496597269717, 24.028657125185283, 23.977118019963747, 23.92014484060037, 23.368047936851802, 23.641608402743394, 22.24382036040257, 21.186770638002955, 22.595994256022106, 21.859099034105807, 21.26534215256648, 21.591354734651826, 21.867605661012256, 21.173676235899062, 22.976958131031367, 21.494660086797374, 21.50084940298445, 22.752545312037594, 20.682160327749465, 21.66009898446791, 22.72954636433276, 21.11959936911907, 21.72106746124736, 22.61729135470916, 21.615026580143002, 21.797275861291954, 20.444952395994513, 22.236018875666367, 20.455047898536094, 20.65298828834646, 20.532172267178627, 22.131252679500847, 21.136604569012537, 21.676002438848656, 21.780911062868217, 20.356429005686156]),
array([23.99896689872301, 21.29074314306712, 15.01623462774724]),
array([23.015735249146516, 21.12246418941679, 21.19650260209131]),
array([23.20484261116118, 23.296635430213215, 23.595120671977835, 21.834483194714267, 22.866744091668057, 22.380670162278022, 21.939482305103333, 22.06893946718261, 22.31746055547692, 22.16631230049523, 22.278595583528453, 22.000517139871512, 22.32147180328683, 22.955648344599982, 22.742791736850215, 22.61983670657284, 21.848501120516172, 22.309492569100648]),
array([22.534253072928024]),
array([21.89162861432136]),
array([22.149686507024068, 22.373476324912133, 22.13400098787761, 22.727994133971205, 22.279389715461633]),
array([21.556842073912982, 21.664521982752614, 22.23379633210613, 21.858922826107765]),
array([22.108079631016505, 22.468354043024632, 22.19697159034776, 22.319471925366418, 22.04987242589144, 22.154149191049772, 21.985809048337327]),
array([21.6746938455199, 21.461352921564334, 21.45635764850666, 22.197340333713086, 21.792525907736895, 21.317649374884166, 20.495431349225928, 21.626706401902197, 20.290638877460438, 19.169755195926765]),
array([21.336080461772255, 21.168325314421985, 21.08617860878633, 20.708363033203224, 21.97307915274237, 21.02334355470829, 20.887786686237217, 20.634403537945346, 21.157665970416, 20.854636848652447, 21.86171909464273, 21.104792935912236, 20.851243443697864, 20.868307007041466, 19.57570579008131, 20.18385953934223, 18.045648116372703, 19.54428212939912, 20.140734854485537]),
array([21.840932150363404, 21.676446149676213]),
array([21.11418262497251, 20.769814354559973, 20.712722576778592, 20.827525687979595, 21.321912612503542]),
array([21.49682175733214, 21.23675315567106]),
array([20.50289092358351, 20.657381122205372, 21.101385004420333, 20.954628090309612, 21.204634193370477]),
array([20.74055189369423, 20.55610733683951, 20.26780942265832]),
array([20.63769290100033, 17.72377829230846, 14.859374291021721, 14.462007100765417, 13.368767238922223, 7.925608230207846, 8.238112433884794, 6.453931283572529, 4.148329887438724, 2.444266199124398, 2.3219964867795673]),
array([20.29242918575033]),
array([18.663560927770877, 19.11258901263367, 17.627458920739777, 16.88326962839121, 17.26445103451868, 17.623947964796372, 14.145908603344093, 13.866707561907658, 14.323327382149696, 15.625103446367325, 14.68145261298908, 14.58677358876799, 15.472426183576665, 14.83387744390769, 14.266516367220065, 15.758926040263233, 14.97324453898306, 13.665587135896907, 13.05379775105201, 13.060836687480393, 12.728326580545758, 13.709665070146565, 12.602632448309588]),
array([14.553427580043614]),
array([19.572407777813083]),
array([18.717478863864415, 17.732686408548485, 19.294321997803245]),
array([18.654407921050733]),
array([18.583515893506817, 18.541988630612497]),
array([14.852483485762587, 15.014155570176895]),
array([16.598642841096726, 17.59963308558488, 17.124960508411423, 17.905269198268222, 17.565673668304736, 16.701450750481612]),
array([13.877123543587231]),
array([16.329359311140006, 16.227827905579147, 15.16907487372206, 13.95823887537551, 15.001823598957449, 15.74899627805131, 15.388303592454418, 14.78704501450538, 15.611608276583805, 15.101831382869353, 15.744275561203166, 15.622218358446483, 13.990263017331165, 14.111494879354556, 14.56960875610538, 14.855969448694708, 14.305600410438325, 14.742865962801954, 14.494779802128583, 15.840761857763306, 15.7073223016059, 13.880979076065058, 14.61165671244071, 15.568596212639207, 15.231980576654976, 14.504360887889897, 15.022645178746405, 13.178019940634478, 12.818405537344894, 13.245584798848885, 11.893440786780241, 12.25200793099402, 12.806391707050054, 12.223812216709721, 11.410069065890742, 8.720632006449396, 9.744420522498427, 9.482342649708244, 9.087281597104418, 8.788587779980398, 5.781451211956919, 6.668935583377433, 6.413142990220242, 4.6215659561532165, 3.9395327088936103, 3.6017593628369173, 5.197230283205585, 3.9941215369646805, 3.8431031542088294, 4.305070042582868, 3.22252376728114, 2.0014370187317243, 2.4549856828418357, 1.8995844570584826, 1.959566614005786, 1.9965056892131186, 2.0348191047017696, 2.5694477383812155, 1.9039847789100333, 2.2611475344697465, 2.4496817160868245, 1.8466715455878444, 1.7449115267894477, 0.8033155088476166, 1.5176342336139885, 1.1717038108070468, 1.115908505167608, 1.0335458805340245, 1.1354319782623468, 1.7616310107077908, 1.598662368408092, 1.39954389818661, 1.5482541955645976, 0.7259902069893227, 0.27519913616146685, 0.3624104812967438, 0.5607990906376373, 0.033108289604995636, 0.0]),
array([14.04852707982705, 13.897318673426915]),
array([15.76149058362325, 14.693180371994302, 15.365372711440923, 14.371436823015419, 15.550733525727862, 14.161470842507025, 15.102851459665523, 15.706240905754882, 15.441866033156938, 15.492814362688627, 15.005470229254094, 13.906542588109449, 14.185593306280493, 15.450916415807814, 13.719416741036364]),
array([14.66757106291379, 15.103763171840322, 12.118576314496927, 7.894898160060533, 6.25466616905201, 4.583852010947178]),
array([16.073020547199835, 16.37743464610627, 16.32886495396262, 15.266415559355362, 15.133930624908675, 15.408899378938239, 15.792578453378335, 15.249605840642907, 14.721572425357884]),
array([16.234217947494162, 15.962782365850334]),
array([15.985710844055596, 15.724132484524134, 13.948134234459967, 14.179672605182915, 13.922332923035174, 15.43484982700789, 15.367428485890391, 15.702892273482833, 15.708495512706588, 15.253633856460063, 15.820301869271871, 14.11867631976557, 15.174600230270256, 14.393733847913206, 14.921528984304793, 15.123854757074724, 15.830273210395376, 14.99804303639564, 15.675570560099427, 14.116158729962484, 14.400648021372717, 15.90372161725638, 14.719261811505225, 13.984220350694455, 15.966704885230552, 13.824830174373124, 13.950607684900046, 15.282463494349713, 13.911310984710937, 12.988929957526533, 12.825860527213527, 12.962801752662987, 13.184447028668885, 12.842471492710823, 11.902051180136834, 13.336305164065008, 13.546614894228364, 13.512364028138078, 12.20961208510454, 13.145143028023183, 12.930515058820923, 11.590106587050856]),
array([16.15109077321979, 15.186183421004882, 14.68121175867325, 15.546764473407768, 13.56641719393987, 13.142556088052796, 7.944312787627291, 6.501205873880508]),
array([15.536588575569317, 15.397039723866571, 15.461393358511327]),
array([15.198593498853707, 14.838560044551514, 14.323133396594239, 14.370658668651906, 15.635016830878328, 12.993928345298837]),
array([15.017342692945592, 7.315691470850965, 10.233968668507064, 10.258194892764603]),
array([14.386865658611532, 14.280118674625884, 13.86405731448085, 11.679739820283318, 11.635088372845669, 13.747922703168094, 13.470415883220943, 8.789805263934378, 10.665806150779806, 9.711859541198057, 7.339880909307881, 10.876481150310415, 8.778941926237842, 10.786295658414268, 10.282776279026294, 6.2082944317243065, 4.430785745344478, 4.911939362717671, 4.278882644153253, 4.0815851669864776, 5.267106389498803, 2.241781745624289, 2.129369835858366, 2.0964921873636793, 2.1595916948185083, 1.6107058616449537, 1.0303386895996702, 1.4731653920422316, 1.0593807262929442]),
array([15.124468495792653, 14.219565708790117]),
array([14.983693340139288, 15.212515099025516, 15.223347472810525, 15.338291530017333]),
array([13.603619989420462, 13.205569273807884, 12.927645850200081, 12.547930726693084]),
array([13.982950925270707, 14.773520540776493, 12.451971312972837, 12.648678338252598]),
array([14.375090368984672, 13.888277252488631, 12.340591586864912, 13.567487977771053, 12.60025473137357, 12.519493890807048, 12.048347871370478, 12.86421596412909, 12.074505860798698, 12.12423142326412, 13.620981032030446, 7.974587597285468, 9.374599973898292, 7.810339590243501, 7.727048097025614, 8.812661462846112, 9.92335723825104, 9.68525041387792, 9.077457356202752, 8.69813398169524, 9.856617855101648, 8.422588771897809, 5.692310130127471, 5.541095724959317, 5.697669361951375, 5.261837580647065, 3.984722538203581, 4.1849360320509765, 4.162749818780114, 4.123010942209256, 4.443830678398433, 4.875110175516097, 2.6804196550160215, 3.0675647050478334, 1.8158600479521443, 1.8588064075142472, 1.9381697283956216, 2.4156945601187867, 1.863949031847005, 2.0248446536358644, 2.1521144320875694, 2.311113060968798, 1.0848901759573146, 0.790780151144175, 1.443640931625853, 1.3359624123478797, 1.2547108631084734, 1.2599917299853276, 1.0071808688562904, 1.537806377079976, 1.4144430393258918, 1.0803545342947456, 1.5457938479438615, 0.907653966047684, 1.384137774765437, 1.5113499196933566, 1.7357151990161408, 0.19015122557998443, 0.34208104886734775, 0.5443365799966697, 0.3283520438082634, 0.7182241479952105, 0.0]),
array([12.591006687767706, 13.051441962162867, 12.370553037147381, 12.570062076989638, 8.860775295010738, 11.039815503943824]),
array([12.834996033201183, 7.532371745766612]),
array([8.07065413860816, 9.20070313303245, 5.294808167707867, 2.42694292337954]),
array([13.404548425785054, 13.446253250779737, 13.422971631825389]),
array([12.093636365871388, 13.479056372369367, 13.429557931240977, 12.782120191753377, 13.471366343371246, 12.378640163687887, 13.061566317540224, 11.93033077708095, 13.467062748614673, 12.463593329693701, 13.151336475835423, 11.680035953289918, 11.526427805620642, 11.33420202041421, 11.21533488934252]),
array([13.22024051456916]),
array([12.916372806835577]),
array([11.103405465816934]),
array([12.130685124263207, 9.92137474993828, 7.925253097069905, 11.38545300819762, 9.193297416825152]),
array([10.214309093222314, 8.073451832175534, 3.4005285200197886, 2.1754817880890283, 2.4048149912689625]),
array([11.74522646344518, 12.050443338324735, 8.351907812750685, 11.3308566131641, 7.6421854293211045, 6.56468032456856, 5.796141542146872, 5.796843907232865, 5.252736913755746]),
array([12.123162378503626, 11.11320168099194, 9.251756935322062]),
array([9.38948632627119, 8.378561191939644, 11.167024551077624, 8.989561953546561, 7.706112721810495, 8.539331910806505, 9.260524253526018, 11.028134836364227, 7.7691534246202565, 5.739007511877678, 5.661266685755695, 7.13980803595757]),
array([10.780817112483572, 10.780225391460995]),
array([9.6667503876043, 8.260988937395949, 11.082998615391542, 10.877122268586474, 6.0053405652986385, 7.075357679575877, 6.785360043962735, 5.145846859671113, 4.9128450157473535, 5.3226658786660925, 3.895060691540094, 4.594292633000376, 3.9370261444372314, 1.8971336922726196, 2.0608855513292896, 1.8779518536467497, 1.8337544281765679, 1.9566272796183708, 1.3704112601467475, 1.6159982398731907, 1.232615515031636, 1.4222166836947312, 0.1737388555206204, 0.0]),
array([10.128574266361074, 7.578131436030276, 8.345502751681982, 6.8603068889319925, 4.401181423003018, 4.5509548587322755, 1.1862570195588098, 0.0]),
array([8.090105290173852, 3.8153627962181673, 1.8376261086968315, 2.1283340871344643, 2.0720980780920306, 2.567432823477319, 1.6690461813028414, 0.0]),
array([10.467799636600763]),
array([9.455185362622286, 8.01698418933615, 9.257842758902164, 7.567649402774736, 7.777943277894023, 8.17028373229023, 9.17828800237552, 8.904679729986913, 8.18231787891808, 9.043929108461862, 8.754422762461074, 7.405600790310045, 6.497962064677299, 6.681930234234792, 5.3892862883765495, 6.879541634416839, 5.6263161165469056, 5.000405324220571, 4.610521177555931, 4.5316273301189725, 3.7451070129217734, 3.945742803915166, 4.543685784557816, 4.669291983836443, 5.298490674662209, 5.13466235862992, 4.689165558269924, 4.0852239323114015, 5.285449913312355, 3.6679033236338667, 4.751147279480151, 2.630472626644742, 2.666892569196217, 2.0834842950338306, 2.253842970528902, 2.3995027987064064, 2.017557088944928, 2.4886889885372625, 2.420525338140266, 1.9892739974534726, 2.111546025366413, 1.984857373081905, 2.20251462554923, 2.462884093528864, 2.033820728640543, 1.9976164232286693, 2.487338966459641, 2.540995319836755, 2.3035203248778924, 2.1122008430110486, 2.3395503458009985, 2.4473603366046133, 1.8795901917765006, 0.7824088107850153, 0.8240326534305006, 0.9380699225203428, 1.3761839713449249, 0.9465504207903039, 1.0796224853817085, 1.218656064730742, 1.0342748731108178, 1.068948220657345, 0.8693666210165273, 0.9331246970202851, 0.9674543214956811, 1.0764332526993203, 1.4356071123986642, 1.3278918157572794, 0.9707915614101436, 1.5955737257433984, 1.7766875421498676, 0.9335042202190938, 1.767336280498664, 0.777878115580828, 0.45618809259716747, 0.4183062455414057, 0.6290709282579412, 0.7014167557935471, 0.7570712408946978, 0.37428225166531837, 0.0]),
array([8.035972074205917, 4.525768832400729, 4.817013322893513, 4.144543472730864, 1.9847109279843922, 0.0]),
array([0.2629848286486457, 0.0]),
array([7.642424501426136, 7.383605867967674, 7.10755003439836, 6.899188543551603, 7.2352787400612515, 6.857391648698834, 6.734248102722328, 6.511533877098948]),
array([9.66448495011544]),
array([8.16457878366702]),
array([8.79360009834965, 8.307182710435743, 4.95179031158931, 4.042028862042091, 2.5130388935892687, 2.132089656009444, 0.0]),
array([8.57546414068112, 8.452693628446806, 7.985684687525647, 8.09193377152569, 5.975080139258384, 6.864167825769181, 6.532599897158472, 6.685941262359492, 5.385730719586002, 7.1672033732233755, 4.793674075026305, 4.55869606901017, 4.455936077624516, 4.98348402428378, 5.0265859786244, 5.016860806969571, 4.360515408385163, 4.324274000545341, 2.3129844039144563, 1.9545050968608635, 2.121080429396607, 1.8654867694059813, 2.055044874348396, 0.8397771901693426, 0.961824153965627, 1.778038421665579, 1.3945002908386879, 1.7968786909510828, 1.4537506209625175, 0.5954951389039453, 0.6308084474531296, 0.2887032493041481, 0.04161490336693621, 0.0]),
array([8.606870691998889]),
array([7.400997256517768, 5.506246105375872, 5.501362324088586, 7.196580079976147, 7.096144392600757, 4.681223202546804, 4.271422357643495, 4.744533459328785, 4.182473651517353, 4.750999292456608, 5.023222090325097, 4.422298555902591, 1.8268847824482477, 2.140335231752084, 2.1152892065998707, 2.179711105400944, 2.519081488467375, 1.4155054191687215, 1.3522743398956587, 1.44761599073395]),
array([8.100703972787846, 8.0155465091615, 7.485518580354659, 7.427201187373287, 7.250292665250215]),
array([6.35531681612983]),
array([7.893281463843355, 4.970469010580868, 4.739422705560411, 4.94951605450371, 4.1817340997015595]),
array([7.865750084451861, 5.743446281263775, 6.232707539674757, 6.175767461710958, 5.962118943346462, 7.101149857822854, 6.771600344736157, 3.8641159066755355, 4.413118800267795, 3.937974845445124, 5.209940728984264, 4.9131706289461485, 3.2491069379908653, 2.108469038764813, 2.490325489400218, 2.3914950972720983, 2.5020454544345716, 2.5295292149421242, 2.5449133328758093, 1.9551045310690098, 1.8772143200223739, 2.5000598962310425, 1.6300700017053502, 1.511412065945729, 1.4627144389202071, 1.5418349501740958]),
array([7.969240680439732, 7.1412426804531925]),
array([6.987471024237083, 3.838992918514979, 4.432717293702133, 0.7736891412239916, 0.0]),
array([7.438619561637877, 7.298649539345064]),
array([6.080547078829998]),
array([4.59266319533875]),
array([6.272037187664282]),
array([6.0741932209615594]),
array([1.39389386040193, 1.3030323786618743, 0.7371730620452585, 0.0]),
array([5.823511940987383, 6.259060295229361]),
array([4.055371506853502, 4.808553288525427, 5.087150290758385, 2.2716434053298746, 2.008864958474445, 2.459285553615064, 2.0107523691759757, 1.6060353520064454, 0.9945413187718056, 1.7205618541328058, 0.83860085261273, 0.0]),
array([5.467298037544903, 4.979621207460034, 3.459672195465503, 2.2534406489137995, 2.4605967024725968, 2.2030476663378344, 1.506950275110403, 1.0917678088460248, 1.0372520152884985, 1.5095607188564681, 0.7936021722579543, 0.23512457058803193, 0.34942333395305236, 0.470506259022605, 0.0]),
array([4.471130476672853, 4.599613162313002, 5.014897864141684]),
array([3.988365438864646, 3.810642420937504, 2.2620817416712984]),
array([5.134857287623208]),
array([4.842613268973255]),
array([3.7175945196009637, 2.096399563727918, 2.455282583033849, 1.0212261410436314, 0.0]),
array([4.6109520184758335, 4.463780824983701, 4.677661034542753, 4.379953452360567, 4.7198597083769265]),
array([4.435002855733484, 3.736721918506742]),
array([2.2414765228374094, 1.9085024620605107, 1.0872710940151835, 1.361307988345425, 0.9618615460080105, 1.2123945436926173, 1.5573286334667125, 0.0]),
array([4.581304190002331, 4.2336162305431575, 4.166480583233454]),
array([3.810529607546166, 3.6340348220649448, 1.8491382012763746, 2.2311723628233704, 1.8942386332776322, 1.857957409381115, 1.0550328547962815, 1.379058312343842, 0.589911500700097, 0.0]),
array([2.0691925412780883, 0.6006093210109444, 0.0]),
array([3.9211312148266466, 3.8991530504982372, 2.0428306701527648, 1.0555559905686247, 1.1869404047276602, 1.301613358652277, 0.9329476195634967, 0.783026107317581, 0.0]),
array([2.582112255586465, 2.4778930345378223, 2.5425590298031584, 1.9273942448096397, 2.1421079165667702, 1.8937693754487304, 2.0150717953943262, 1.9213693755858232, 2.0903615515880754, 1.0642119672099426, 1.5370161002841962, 1.2522922960166993, 1.7849427344355673, 1.145574221122288, 1.1712886452982845, 0.5743923187416697, 0.7153370803717615, 0.0]),
array([2.1853573238388817, 0.1344001071618829, 0.0]),
array([3.050625087185398, 3.439634351864053, 2.3349819324769934, 1.9390850671811592, 2.1654974882330738, 2.231931138716919, 1.9312964722309136, 2.33052359281348, 2.244139958013123, 2.0833259673082405, 2.3224696481189797, 1.2682200882662142, 1.6603868399896866, 1.2274860276964932, 1.1544231601522346, 1.3668308930944708, 0.7940856645748873, 1.19878875917103, 1.7575774759175218, 0.15776565052083857, 0.35670110162994934, 0.6605691758558175]),
array([2.359158933259809, 0.0]),
array([2.036958624593759, 1.9210324000083565, 1.7078040356549318, 1.5788392094569121, 0.0]),
array([2.344184607813156, 0.0]),
array([2.809447240356695, 2.7164190705804083, 2.7817070881023453, 2.1097950658771554, 2.3065260126448672, 2.2208613795850742, 2.1953186352854352, 2.3277279602239176, 2.1929553243378135, 1.3991491189082577, 0.9669620038599144, 1.2543676382318478, 1.1695474316519119, 1.131651361825237, 1.581170704896065, 1.0574898549234253, 1.666560796096396, 1.691727267084957, 1.792611475507274, 0.7137903859638026]),
array([2.946930998137071, 2.02014661801536, 2.514064357703678, 2.540638450417061, 2.104929452107929, 2.4434706985734094, 2.0006235858020904, 2.2117156362432913, 1.7725552027151896]),
array([2.281849507526382, 2.0511400280174894, 1.915646916379244, 1.0304185656892788]),
array([2.7120278875086705, 2.2518892857290718, 2.558100865239153, 2.4883671634994156, 2.00159448943836, 1.7293134260451903, 1.6408355478775367, 1.1337645357710673, 1.3416584107368772, 1.0154075074788365, 0.5515042606661019, 0.42562406900412114, 0.6104353776944633, 0.12027083394576928, 0.0]),
array([2.7965789515386152, 2.3086764106682605, 2.0606718716927044, 2.5089999455642236, 1.276644249251875, 1.7203293477539605, 0.8448622423769205, 0.8390097239974602, 0.8563658608672465, 1.6110604483991577, 1.1391541642777283, 1.624008599650527, 0.47217137039216733, 0.30299091678101187, 0.0]),
array([2.1089565774427546, 2.0150518021131236, 2.3141874149024013, 2.260759757721158, 1.812542282252357, 1.1643150457452567, 1.277831016968813, 1.0871496604850162, 0.9005124295316309, 1.242187929498599, 1.5690148527248953, 0.23543660287851143, 0.0]),
array([1.8228851213036532, 1.83516288829765, 2.2306349679751953, 2.3299189017711592, 2.0090045848391607, 2.1623047185075652, 2.2816119574606, 2.201954228974451, 2.481392995498542, 2.2962601667712845, 2.3973100263867253, 2.084601957872425, 2.2680135698148263, 2.267939756915744, 2.5694300607415803, 1.6840001967549214, 1.035656169360854, 1.566704363571781, 1.3141732669918729, 1.387568256308304, 1.388945852145714, 1.2405593791892322, 1.6198758643305424, 1.6443134720497319, 0.9965602072097299, 0.8654141640242711, 0.21291541543356185, 0.650008980897503, 0.7657838253882395, 0.700672906021664, 0.557363095389237, 0.0]),
array([2.326430600296911]),
array([2.020766192750363, 1.868003966226719, 2.2781132513195264, 0.8063895412929781, 0.8230290015682001, 0.2693560312709513, 0.0]),
array([2.2139302988940917, 1.8667164687242501, 1.8983949493406682, 1.8565302549011637, 0.9938881096474109, 1.0673094897387831, 0.8657777700452246, 1.2904147959940035, 0.5538594393697875, 0.4607818747435097, 0.0]),
array([1.1023397192172986, 1.1350528749496172, 1.6986753776551708, 0.0]),
array([0.7483748841631604, 0.0]),
array([1.7856297642979968, 1.5701282047789615, 0.8383675884956506, 0.5674471597089492, 0.289511788426921, 0.496204766051226, 0.0]),
array([0.8611575859366131, 1.5305079046472985, 0.9978603341441187, 1.1663045922579336, 1.3861613771175287, 1.0796635530941827, 1.5358308775937495, 1.5540373643541106]),
array([1.4191695031649179, 0.9944998460247236, 1.185879870956128, 0.9589765987784135, 0.8426419685929214, 1.1769978476623066, 1.4016117023543568, 0.9161705875596845, 1.1626368295791216, 0.9750993764177903, 0.9980290937551062, 1.1505414836668149, 0.28477793422826475, 0.5775948373583806, 0.0]),
array([1.4625369438774862, 1.072747138886978, 0.0]),
array([1.3235173595272043, 1.1246777512640398, 1.0275549936917703, 1.2017546676906568, 1.021032020817193, 0.783076367523712, 0.9669830638968626, 0.2861204165404577, 0.6642941208811044, 0.0]),
array([0.22847929284380175, 0.0]),
array([0.5815253184436889, 0.0]),
array([0.839874135752239, 1.2393787993525496, 1.1046755779016149, 0.7285273700490633, 0.0]),
array([0.9065832572193055, 0.9503172990997298, 0.8783615045158732, 0.5289228675766245, 0.18592180008832804, 0.3249955837097089, 0.152208323517113, 0.0]),
array([0.9413443114889095, 0.8800175459441111, 0.4295472581650248, 0.02906260703997303, 0.0]),
array([0.5326168553816513, 0.0]),
array([0.18309448980260934, 0.0]),
array([0.18312308167620883, 0.5992947684722558, 0.0]),
array([0.4446522290346282, 0.31409153871320705, 0.0]),
array([0.3681778889433391, 0.2800763561099843, 0.0]),
array([0.39767925709840446, 0.0]),
array([0.2915905582787194, 0.2620310360353328, 0.21996144717856486, 0.0]),
array([0.1737846210142348]),
array([0.0817040413580766, 0.0]),
array([0.026772080347671465, 0.0])
]
d = [data_1]
names = ["36"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T1', 'T2', 'T4', 'T5', 'T6', 'T8', 'T9', 'T10', 'T12', 'T13', 'T14', 'T15', 'T16', 'T18', 'T19', 'T20', 'T22', 'T23', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T42', 'T44', 'T45', 'T46', 'T47', 'T50', 'T51', 'T52', 'T53', 'T54', 'T56', 'T57', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T66', 'T68', 'T69', 'T70', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T80', 'T82', 'T83', 'T84', 'T87', 'T89', 'T90', 'T91', 'T94', 'T97', 'T98', 'T99', 'T104', 'T106', 'T107', 'T109', 'T110', 'T111', 'T113', 'T114', 'T115', 'T116', 'T117', 'T119', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T130', 'T131', 'T134', 'T135', 'T137', 'T139', 'T140', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T158', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T168', 'T169', 'T171', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'T179', 'T180', 'T181', 'T183', 'T184', 'T185', 'T186', 'T187', 'T188', 'T189', 'T191', 'T192', 'T193', 'T194', 'T196', 'T197', 'T199', 'T200', 'T203', 'T205']
def get_taxa_names(): return taxa_names