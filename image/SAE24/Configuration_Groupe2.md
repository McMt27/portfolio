
# sw3
!
version 15.0
no service pad
service timestamps debug datetime msec
service timestamps log datetime msec
no service password-encryption
!
hostname SW3-Server
!
boot-start-marker
boot-end-marker
!
!
username admin privilege 15 secret 5 $1$s2Oi$3ZEOKgMQupkhDmWo/THiL0
no aaa new-model
system mtu routing 1500
!
!
no ip domain-lookup
ip domain-name pellet.com
!
!         
crypto pki trustpoint TP-self-signed-3608219008
 enrollment selfsigned
 subject-name cn=IOS-Self-Signed-Certificate-3608219008
 revocation-check none
 rsakeypair TP-self-signed-3608219008
!
!
crypto pki certificate chain TP-self-signed-3608219008
!
!
!
!
!
!
spanning-tree mode mst
spanning-tree extend system-id
!
spanning-tree mst configuration
 name REGION-PELLET
 revision 1
 instance 1 vlan 600, 650, 700
 instance 2 vlan 120, 130, 140, 150, 160
!         
!
vlan internal allocation policy ascending
!
ip ssh authentication-retries 4
ip ssh version 2
!
!
!
!
!
interface FastEthernet0/1
 description PORT_TEST
 switchport access vlan 600
 switchport mode access
!
interface FastEthernet0/2
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/3
 description NON_UTILISE
 shutdown
!         
interface FastEthernet0/4
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/5
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/6
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/7
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/8
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/9
 description NON_UTILISE
 shutdown 
!
interface FastEthernet0/10
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/11
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/12
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/13
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/14
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/15
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/16
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/17
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/18
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/19
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/20
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/21
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/22
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/23
 description NON_UTILISE
 shutdown
!
interface FastEthernet0/24
 description VERS_PROXMOX
 switchport access vlan 700
 switchport trunk native vlan 999
 switchport mode access
 spanning-tree portfast
!
interface GigabitEthernet0/1
 description VERS_SWR1
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
!
interface GigabitEthernet0/2
 description VERS_SWR2
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
!
interface Vlan1
 no ip address
!
interface Vlan600
 description Management_SW3
 ip address 172.16.0.10 255.255.255.0
!
interface Vlan700
 description VERS_PROXMOX
 no ip address
!
ip default-gateway 172.16.0.254
ip http server
ip http secure-server
access-list 10 permit 172.16.0.0 0.0.0.255
!         
vstack
!
line con 0
line vty 0 4
 access-class 10 in
 exec-timeout 5 0
 login local
 transport input ssh
line vty 5 15
 access-class 10 in
 exec-timeout 5 0
 login local
 transport input ssh
!
end



# SWR1
!
version 12.2
no service pad
service timestamps debug datetime msec
service timestamps log datetime msec
no service password-encryption
!
hostname SWR1
!
boot-start-marker
boot-end-marker
!
!
username admin privilege 15 secret 5 $1$4RcJ$192EeKPottrrlwgqAiBvr/
!
!
no aaa new-model
switch 1 provision ws-c3750g-24t
switch 2 provision ws-c3750g-24t
system mtu routing 1500
ip routing
ip domain-name pellet.com
ip dhcp excluded-address 172.16.0.1 172.16.0.20
ip dhcp excluded-address 172.16.0.252 172.16.0.254
ip dhcp excluded-address 172.16.0.1 172.16.0.30
ip dhcp excluded-address 172.16.0.250 172.16.0.254
!
ip dhcp pool POOL-ADMIN-INFO
   network 172.16.0.0 255.255.255.0
   default-router 172.16.0.254 
   dns-server 172.16.2.7 
!
ip dhcp pool POOL-COMPTA-VLAN-120
   network 172.16.3.0 255.255.255.192
   default-router 172.16.3.62 
   dns-server 172.16.2.7 
!
ip dhcp pool POOL-TERRAIN-VLAN-130
   network 172.16.3.128 255.255.255.128
   default-router 172.16.3.254 
   dns-server 172.16.2.7 
!
ip dhcp pool POOL-ADMINISTRATIF-VLAN-140
   network 172.16.3.64 255.255.255.192
   default-router 172.16.3.126 
   dns-server 172.16.2.7 
!
ip dhcp pool POOL-PROD-VLAN-150
   network 172.16.4.0 255.255.255.128
   default-router 172.16.4.126 
   dns-server 172.16.2.7 
!
ip dhcp pool POOL-LOGIS-VLAN-160
   network 172.16.4.128 255.255.255.192
   default-router 172.16.4.190 
   dns-server 172.16.2.7 
!
!
!
!
crypto pki trustpoint TP-self-signed-3375198592
 enrollment selfsigned
 subject-name cn=IOS-Self-Signed-Certificate-3375198592
 revocation-check none
 rsakeypair TP-self-signed-3375198592
!
!
!         
!
!
!
spanning-tree mode mst
spanning-tree extend system-id
!
spanning-tree mst configuration
 name REGION-PELLET
 revision 1
 instance 1 vlan 600, 650, 700
 instance 2 vlan 120, 130, 140, 150, 160
!
spanning-tree mst 0-1 priority 4096
spanning-tree mst 2 priority 8192
!
vlan internal allocation policy ascending
!
ip ssh time-out 60
ip ssh authentication-retries 4
ip ssh version 2
!
!
!         
interface Port-channel1
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
!
interface Port-channel2
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
!
interface Port-channel10
 description TRUNK_INTER_COEUR
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
!
interface GigabitEthernet1/0/1
 description VERS_R1_G0/1
 no switchport
 ip address 172.16.5.2 255.255.255.252
!
interface GigabitEthernet1/0/2
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
 shutdown
 channel-protocol lacp
!
interface GigabitEthernet1/0/3
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
 shutdown
 channel-protocol lacp
!
interface GigabitEthernet1/0/4
!
interface GigabitEthernet1/0/5
 description VERS_SW3
!
interface GigabitEthernet1/0/6
!
interface GigabitEthernet1/0/7
!
interface GigabitEthernet1/0/8
!
interface GigabitEthernet1/0/9
!
interface GigabitEthernet1/0/10
 description TEST_COMPTA
 switchport access vlan 120
 switchport mode access
 spanning-tree portfast
!
interface GigabitEthernet1/0/11
 description TEST_TERRAIN
 switchport access vlan 130
 switchport mode access
 spanning-tree portfast
!
interface GigabitEthernet1/0/12
 description TEST_ADMIN
 switchport access vlan 140
 switchport mode access
 spanning-tree portfast
!         
interface GigabitEthernet1/0/13
 description TEST_PROD
 switchport access vlan 150
 switchport mode access
 spanning-tree portfast
!
interface GigabitEthernet1/0/14
 description TEST_ADMIN_INFO
 switchport access vlan 600
 switchport mode access
 spanning-tree portfast
!
interface GigabitEthernet1/0/15
!
interface GigabitEthernet1/0/16
!
interface GigabitEthernet1/0/17
!
interface GigabitEthernet1/0/18
!
interface GigabitEthernet1/0/19
!
interface GigabitEthernet1/0/20
!
interface GigabitEthernet1/0/21
 channel-group 1 mode active
!
interface GigabitEthernet1/0/22
 channel-group 1 mode active
!
interface GigabitEthernet1/0/23
 channel-group 2 mode active
!
interface GigabitEthernet1/0/24
 channel-group 2 mode active
!
interface GigabitEthernet2/0/1
 description UPLINK_SW1
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
 channel-protocol lacp
 channel-group 1 mode active
!
interface GigabitEthernet2/0/2
 description UPLINK_SW1
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
 channel-protocol lacp
 channel-group 1 mode active
!
interface GigabitEthernet2/0/3
 description UPLINK_SW2
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
 channel-protocol lacp
 channel-group 2 mode active
!
interface GigabitEthernet2/0/4
 description UPLINK_SW2
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
 channel-protocol lacp
 channel-group 2 mode active
!
interface GigabitEthernet2/0/5
 description VERS_SW3_GI0/1
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
!
interface GigabitEthernet2/0/6
!
interface GigabitEthernet2/0/7
 description LIEN_ISL_VERS_SWR2
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
 channel-protocol lacp
 channel-group 10 mode active
!
interface GigabitEthernet2/0/8
 description LIEN_ISL_VERS_SWR2
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
 channel-protocol lacp
 channel-group 10 mode active
!
interface GigabitEthernet2/0/9
!
interface GigabitEthernet2/0/10
!
interface GigabitEthernet2/0/11
!
interface GigabitEthernet2/0/12
!
interface GigabitEthernet2/0/13
!
interface GigabitEthernet2/0/14
!
interface GigabitEthernet2/0/15
!
interface GigabitEthernet2/0/16
!
interface GigabitEthernet2/0/17
!
interface GigabitEthernet2/0/18
!
interface GigabitEthernet2/0/19
!
interface GigabitEthernet2/0/20
!
interface GigabitEthernet2/0/21
!
interface GigabitEthernet2/0/22
!
interface GigabitEthernet2/0/23
!
interface GigabitEthernet2/0/24
!
interface Vlan1
 no ip address
 shutdown
!
interface Vlan120
 ip address 172.16.3.60 255.255.255.192
 ip access-group ACL-METIERS in
 standby 120 ip 172.16.3.62
 standby 120 priority 110
 standby 120 preempt
!
interface Vlan130
 ip address 172.16.3.252 255.255.255.128
 ip access-group ACL-METIERS in
 standby 130 ip 172.16.3.254
 standby 130 priority 110
 standby 130 preempt
!
interface Vlan140
 ip address 172.16.3.124 255.255.255.192
 standby 140 ip 172.16.3.126
 standby 140 priority 110
 standby 140 preempt
!
interface Vlan150
 ip address 172.16.4.124 255.255.255.128
 ip access-group ACL-METIERS in
 standby 150 ip 172.16.4.126
 standby 150 priority 110
 standby 150 preempt
!
interface Vlan160
 ip address 172.16.4.188 255.255.255.192
 ip access-group ACL-METIERS in
 standby 160 ip 172.16.4.190
 standby 160 priority 110
 standby 160 preempt
!
interface Vlan600
 ip address 172.16.0.252 255.255.255.0
 ip access-group ACL-ADMIN-INFO in
 standby 60 ip 172.16.0.254
 standby 60 priority 110
 standby 60 preempt
!
interface Vlan700
 ip address 172.16.2.124 255.255.255.128
 standby 70 ip 172.16.2.126
 standby 70 priority 110
 standby 70 preempt
!
ip classless
ip route 0.0.0.0 0.0.0.0 172.16.5.1
ip http server
ip http secure-server
!
!
ip access-list extended ACL-ADMIN-INFO
 permit tcp any 172.16.0.0 0.0.255.255 eq 22
 permit ip any 172.16.2.0 0.0.0.127
 deny   ip any 172.16.0.0 0.0.255.255
 permit ip any any
ip access-list extended ACL-METIERS
 permit ip any 172.16.2.0 0.0.0.127
 permit ip any 172.16.3.64 0.0.0.63
 deny   ip any 172.16.0.0 0.0.255.255
 permit ip any any
!
access-list 10 permit 172.16.0.0 0.0.0.255
!
!
!
line con 0
line vty 0 4
 access-class 10 in
 exec-timeout 5 0
 login local
 transport input ssh
line vty 5 15
 access-class 10 in
 exec-timeout 5 0
 login local
 transport input ssh
!
end



# SWR2
!
version 12.2
no service pad
service timestamps debug datetime msec
service timestamps log datetime msec
no service password-encryption
!
hostname SWR2
!
boot-start-marker
boot-end-marker
!
!
username admin privilege 15 secret 5 $1$F9wP$bkUDTZEvqFA3JmRKVFN3b1
!
!
no aaa new-model
switch 1 provision ws-c3750g-16td
switch 2 provision ws-c3750g-24ts-1u
system mtu routing 1500
ip routing
ip domain-name pellet.com
!
!
!
!
crypto pki trustpoint TP-self-signed-2461536384
 enrollment selfsigned
 subject-name cn=IOS-Self-Signed-Certificate-2461536384
 revocation-check none
 rsakeypair TP-self-signed-2461536384
!
!
!
!
!
!
spanning-tree mode mst
spanning-tree extend system-id
!
spanning-tree mst configuration
 name REGION-PELLET
 revision 1
 instance 1 vlan 600, 650, 700
 instance 2 vlan 120, 130, 140, 150, 160
!
spanning-tree mst 1 priority 8192
spanning-tree mst 2 priority 4096
!
vlan internal allocation policy ascending
!
ip ssh time-out 60
ip ssh version 2
!
!
!
interface Port-channel3
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
!
interface Port-channel4
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
!
interface Port-channel10
 description TRUNK_INTER_COEUR
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
!
interface GigabitEthernet1/0/1
 description VERS_R1_G0/2
 no switchport
 ip address 172.16.5.6 255.255.255.252
!
interface GigabitEthernet1/0/2
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
 shutdown
 channel-protocol lacp
!
interface GigabitEthernet1/0/3
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
 shutdown
 channel-protocol lacp
!
interface GigabitEthernet1/0/4
!
interface GigabitEthernet1/0/5
!
interface GigabitEthernet1/0/6
!
interface GigabitEthernet1/0/7
!
interface GigabitEthernet1/0/8
!
interface GigabitEthernet1/0/9
!
interface GigabitEthernet1/0/10
!
interface GigabitEthernet1/0/11
!
interface GigabitEthernet1/0/12
!
interface GigabitEthernet1/0/13
!         
interface GigabitEthernet1/0/14
!
interface GigabitEthernet1/0/15
!
interface GigabitEthernet1/0/16
!
interface TenGigabitEthernet1/0/1
!
interface GigabitEthernet2/0/1
 description UPLINK_SW1
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
 channel-protocol lacp
 channel-group 3 mode active
!
interface GigabitEthernet2/0/2
 description UPLINK_SW1
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
 channel-protocol lacp
 channel-group 3 mode active
!
interface GigabitEthernet2/0/3
 description UPLINK_SW2
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
 channel-protocol lacp
 channel-group 4 mode active
!
interface GigabitEthernet2/0/4
 description UPLINK_SW2
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
 channel-protocol lacp
 channel-group 4 mode active
!
interface GigabitEthernet2/0/5
 description VERS_SW3_GI0/2
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport mode trunk
!         
interface GigabitEthernet2/0/6
!
interface GigabitEthernet2/0/7
 description LIEN_ISL_VERS_SWR1
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
 channel-protocol lacp
 channel-group 10 mode active
!
interface GigabitEthernet2/0/8
 description LIEN_ISL_VERS_SWR1
 switchport trunk encapsulation dot1q
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120,130,140,150,160,600,650,700,999
 switchport mode trunk
 channel-protocol lacp
 channel-group 10 mode active
!
interface GigabitEthernet2/0/9
!
interface GigabitEthernet2/0/10
!
interface GigabitEthernet2/0/11
!
interface GigabitEthernet2/0/12
!
interface GigabitEthernet2/0/13
!
interface GigabitEthernet2/0/14
!
interface GigabitEthernet2/0/15
!
interface GigabitEthernet2/0/16
!
interface GigabitEthernet2/0/17
!
interface GigabitEthernet2/0/18
!
interface GigabitEthernet2/0/19
!
interface GigabitEthernet2/0/20
!
interface GigabitEthernet2/0/21
!         
interface GigabitEthernet2/0/22
!
interface GigabitEthernet2/0/23
!
interface GigabitEthernet2/0/24
!
interface GigabitEthernet2/0/25
!
interface GigabitEthernet2/0/26
!
interface GigabitEthernet2/0/27
!
interface GigabitEthernet2/0/28
!
interface Vlan1
 no ip address
 shutdown
!
interface Vlan120
 description Comptabilite
 ip address 172.16.3.61 255.255.255.192
 ip access-group ACL-METIERS in
 standby 120 ip 172.16.3.62
 standby 120 priority 90
 standby 120 preempt
!
interface Vlan130
 description Equipe_Terrain
 ip address 172.16.3.253 255.255.255.128
 ip access-group ACL-METIERS in
 standby 130 ip 172.16.3.254
 standby 130 priority 90
 standby 130 preempt
!
interface Vlan140
 description Administratif
 ip address 172.16.3.125 255.255.255.192
 standby 140 ip 172.16.3.126
 standby 140 priority 90
 standby 140 preempt
!
interface Vlan150
 description Production
 ip address 172.16.4.125 255.255.255.128
 ip access-group ACL-METIERS in
 standby 150 ip 172.16.4.126
 standby 150 priority 90
 standby 150 preempt
!
interface Vlan160
 description Logistique
 ip address 172.16.4.189 255.255.255.192
 ip access-group ACL-METIERS in
 standby 160 ip 172.16.4.190
 standby 160 priority 90
 standby 160 preempt
!
interface Vlan600
 description Admin_Info
 ip address 172.16.0.253 255.255.255.0
 ip access-group ACL-ADMIN-INFO in
 standby 60 ip 172.16.0.254
 standby 60 priority 90
 standby 60 preempt
!
interface Vlan700
 description Serveurs
 ip address 172.16.2.125 255.255.255.128
 standby 70 ip 172.16.2.126
 standby 70 priority 90
 standby 70 preempt
!
ip classless
ip route 0.0.0.0 0.0.0.0 172.16.5.5
ip http server
ip http secure-server
!
!
ip access-list extended ACL-ADMIN-INFO
 permit tcp any 172.16.0.0 0.0.255.255 eq 22
 permit ip any 172.16.2.0 0.0.0.127
 deny   ip any 172.16.0.0 0.0.255.255
 permit ip any any
ip access-list extended ACL-METIERS
 permit ip any 172.16.2.0 0.0.0.127
 permit ip any 172.16.3.64 0.0.0.63
 deny   ip any 172.16.0.0 0.0.255.255
 permit ip any any
!
access-list 10 permit 172.16.0.0 0.0.0.255
!
!         
!
line con 0
line vty 0 4
 access-class 10 in
 exec-timeout 5 0
 login local
 transport input ssh
line vty 5 15
 access-class 10 in
 exec-timeout 5 0
 login local
 transport input ssh
!
end


# R1
!
version 15.7
service timestamps debug datetime msec
service timestamps log datetime msec
no service password-encryption
!
hostname R1
!
boot-start-marker
boot-end-marker
!
!
!
no aaa new-model
memory-size iomem 15
!
!
!
!
!
!
!         
!
!
!
!
!
!
!
!
!


!
!
!
!
no ip domain lookup
ip domain name pellet.com
ip cef
no ipv6 cef
!
multilink bundle-name authenticated
!
!         
!
!
!
!
!
!
voice-card 0
!
!
!
!
!
!
!
!
vxml logging-tag
license udi pid CISCO2911/K9 sn FCZ163460SS
license accept end user agreement
license boot module c2900 technology-package securityk9
license boot module c2900 technology-package uck9
license boot module c2900 technology-package datak9
!
!         
username admin privilege 15 secret 5 $1$EPuU$StnupBRh51svO7SMtr5cv/
!
redundancy
!
!
!
!
!
! 
!
!
!
!
!
!
!
!
!
interface Embedded-Service-Engine0/0
 no ip address
 shutdown
!
interface GigabitEthernet0/0
 description VERS_BOX_ADSL
 ip address 192.168.222.254 255.255.255.0
 ip nat outside
 ip virtual-reassembly in
 duplex auto
 speed auto
!
interface GigabitEthernet0/1
 description VERS_SWR1_TRANSIT
 ip address 172.16.5.1 255.255.255.252
 ip nat inside
 ip virtual-reassembly in
 duplex auto
 speed auto
!
interface GigabitEthernet0/2
 description VERS_SWR2_TRANSIT
 ip address 172.16.5.5 255.255.255.252
 ip nat inside
 ip virtual-reassembly in
 duplex auto
 speed auto
!         
interface Serial0/0/0
 no ip address
 shutdown
 clock rate 2000000
!
interface Serial0/0/1
 no ip address
 shutdown
 clock rate 2000000
!
!
ip forward-protocol nd
!
no ip http server
no ip http secure-server
!
ip nat inside source list 1 interface GigabitEthernet0/0 overload
ip route 0.0.0.0 0.0.0.0 192.168.222.1
ip route 172.16.0.0 255.255.248.0 172.16.5.2
ip route 172.16.0.0 255.255.248.0 172.16.5.6 10
ip ssh version 2
!
ipv6 ioam timestamp
!
!
access-list 1 permit 172.16.0.0 0.0.7.255
access-list 10 permit 172.16.0.0 0.0.0.255
access-list 10 permit 172.16.5.0 0.0.0.7
!
control-plane
!
 !
 !
 !
 !
!
mgcp behavior rsip-range tgcp-only
mgcp behavior comedia-role none
mgcp behavior comedia-check-media-src disable
mgcp behavior comedia-sdp-force disable
!
mgcp profile default
!
!
!
!         
!
!
!
gatekeeper
 shutdown
!
!
line con 0
line aux 0
line 2
 no activation-character
 no exec
 transport preferred none
 transport output lat pad telnet rlogin lapb-ta mop udptn v120 ssh
 stopbits 1
line vty 0 4
 access-class 10 in
 login local
 transport input ssh
!
scheduler allocate 20000 1000
!
end     
# SW1
!
version 15.0
no service pad
service timestamps debug datetime msec
service timestamps log datetime msec
service password-encryption
!
hostname SW1
!
boot-start-marker
boot-end-marker
!
!
username admin privilege 15 secret 5 $1$akER$ojo.KWqIqqkG048FWCmAJ.
no aaa new-model
system mtu routing 1500
!
!
no ip domain-lookup
ip domain-name pellet.com
!
!         
crypto pki trustpoint TP-self-signed-3195865984
 enrollment selfsigned
 subject-name cn=IOS-Self-Signed-Certificate-3195865984
 revocation-check none
 rsakeypair TP-self-signed-3195865984
!
!
crypto pki certificate chain TP-self-signed-3195865984
 certificate self-signed 01
  3082022B 30820194 A0030201 02020101 300D0609 2A864886 F70D0101 05050030 
  31312F30 2D060355 04031326 494F532D 53656C66 2D536967 6E65642D 43657274 
  69666963 6174652D 33313935 38363539 3834301E 170D3933 30333031 30303030 
  35395A17 0D323030 31303130 30303030 305A3031 312F302D 06035504 03132649 
  4F532D53 656C662D 5369676E 65642D43 65727469 66696361 74652D33 31393538 
  36353938 3430819F 300D0609 2A864886 F70D0101 01050003 818D0030 81890281 
  8100A10F 7CDA6E41 69C425D2 3680CE8A 3C689349 9E1C2A43 8B65CF65 1F8605AB 
  B6D430F6 D6D89731 BB35534F 20FFFE44 34E98602 2096DD72 9EA80B6E 07BB0EF2 
  8644AB7F 45EEBB6A F339CEF3 C1D42E28 F62E86FB 17A1E1A7 BD1A0A5C 913449C9 
  9E5EE8ED 4D07C001 AFD67EA3 62FBA851 DD191BF6 DAD6AB6A AD2A41A5 EF4BA200 
  C4F90203 010001A3 53305130 0F060355 1D130101 FF040530 030101FF 301F0603 
  551D2304 18301680 14BD61CD 4013714F E50DB1A4 80DCEEEF 3571EA2F 91301D06 
  03551D0E 04160414 BD61CD40 13714FE5 0DB1A480 DCEEEF35 71EA2F91 300D0609 
  2A864886 F70D0101 05050003 8181002A C38D11C5 F141F44E 5054CEFB BC4776E0 
  C1C3B91F A86A20D7 1AF07DFA 690BBE13 D751F87B 498861B8 216EA023 85DB9C3A 
  5D1D315F 6A534064 AB28AFE5 46BDC8BC 8B588216 798050E8 2C86D871 FE248355 
  42A8E6BD 24B383E1 BBE12BE5 11D9505A D9EFCE78 1B8DEDB2 A0D4BB67 CB7677BE 
  B1EB5135 BE357508 36CF1650 D3B972
  	quit
!
!
!
!
!
!
spanning-tree mode mst
spanning-tree extend system-id
!
spanning-tree mst configuration
 name REGION-PELLET
 revision 1
 instance 1 vlan 600, 650, 700
 instance 2 vlan 120, 130, 140, 150, 160
!
!
vlan internal allocation policy ascending
!         
ip ssh time-out 60
ip ssh version 2
!
!
!
!
!
interface Port-channel1
 description Vers_SWR1_Principal
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
!
interface Port-channel2
 description Vers_SWR2_Secours
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
!
interface FastEthernet0/1
 description TEST_COMPTA_V120
 switchport access vlan 120
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/2
 description TEST_COMPTA_V120
 switchport access vlan 120
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/3
 description TEST_COMPTA_V120
 switchport access vlan 120
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/4
 description TEST_COMPTA_V120
 switchport access vlan 120
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/5
 description TEST_TERRAIN_V130
 switchport access vlan 130
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/6
 description TEST_TERRAIN_V130
 switchport access vlan 130
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/7
 description TEST_TERRAIN_V130
 switchport access vlan 130
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/8
 description TEST_TERRAIN_V130
 switchport access vlan 130
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/9
 description TEST_ADMIN_V140
 switchport access vlan 140
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/10
 description TEST_ADMIN_V140
 switchport access vlan 140
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/11
 description TEST_ADMIN_V140
 switchport access vlan 140
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/12
 description TEST_ADMIN_V140
 switchport access vlan 140
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/13
 description TEST_PROD_V150
 switchport access vlan 150
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/14
 description TEST_PROD_V150
 switchport access vlan 150
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/15
 description TEST_PROD_V150
 switchport access vlan 150
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/16
 description TEST_PROD_V150
 switchport access vlan 150
 switchport mode access
 spanning-tree portfast
!         
interface FastEthernet0/17
 description TEST_LOGIS_V160
 switchport access vlan 160
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/18
 description TEST_LOGIS_V160
 switchport access vlan 160
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/19
 description TEST_ADMIN_INFO_V600
 switchport access vlan 600
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/20
!
interface FastEthernet0/21
 description Uplink_SWR2
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
 channel-protocol lacp
 channel-group 2 mode active
!
interface FastEthernet0/22
 description Uplink_SWR2
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
 channel-protocol lacp
 channel-group 2 mode active
!
interface FastEthernet0/23
 description Uplink_SWR1
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
 channel-protocol lacp
 channel-group 1 mode active
!
interface FastEthernet0/24
 description Uplink_SWR1
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
 channel-protocol lacp
 channel-group 1 mode active
!
interface GigabitEthernet0/1
 description LIEN_VERS_COEUR_SWR
!
interface GigabitEthernet0/2
!
interface Vlan1
 no ip address
 shutdown
!
interface Vlan600
 description Management_SW1
 ip address 172.16.0.12 255.255.255.0
!
ip default-gateway 172.16.0.254
ip http server
ip http secure-server
access-list 1 permit 172.16.0.0 0.0.7.255
access-list 10 permit 172.16.0.0 0.0.0.255
!
vstack
!
line con 0
line vty 0 4
 access-class 10 in
 exec-timeout 5 0
 login local
 transport input ssh
line vty 5 15
 access-class 10 in
 exec-timeout 5 0
 login local
 transport input ssh
!
end


# SW2
!
version 15.0
no service pad
service timestamps debug datetime msec
service timestamps log datetime msec
service password-encryption
!
hostname SW2
!
boot-start-marker
boot-end-marker
!
!
username admin privilege 15 secret 5 $1$97fL$FPeqbuq20LN6x5CchFhM10
no aaa new-model
system mtu routing 1500
!
!
no ip domain-lookup
ip domain-name pellet.com
!
!         
crypto pki trustpoint TP-self-signed-1458781440
 enrollment selfsigned
 subject-name cn=IOS-Self-Signed-Certificate-1458781440
 revocation-check none
 rsakeypair TP-self-signed-1458781440
!
!
crypto pki certificate chain TP-self-signed-1458781440
 certificate self-signed 01
  3082022B 30820194 A0030201 02020101 300D0609 2A864886 F70D0101 05050030 
  31312F30 2D060355 04031326 494F532D 53656C66 2D536967 6E65642D 43657274 
  69666963 6174652D 31343538 37383134 3430301E 170D3933 30333031 30303031 
  30315A17 0D323030 31303130 30303030 305A3031 312F302D 06035504 03132649 
  4F532D53 656C662D 5369676E 65642D43 65727469 66696361 74652D31 34353837 
  38313434 3030819F 300D0609 2A864886 F70D0101 01050003 818D0030 81890281 
  8100CE62 F10BA58F A1266A67 708AF66D 0A55922F 21BC77A0 5C45BF8A 1E1DEE61 
  AB8A41DD C8C64790 C2F28EDF 718EBD43 00BB346B 04119689 86C8FCA6 3A193D3D 
  DDBCEFB0 A62FC417 068C6812 1CCE5EEB 3AEAEA3C 963EAD03 AB1EBFD1 E066333E 
  0369771D 18084721 2FAD9300 8343AF46 CEB48E2D 762CF15A 63DC1272 EC43DC4D 
  53650203 010001A3 53305130 0F060355 1D130101 FF040530 030101FF 301F0603 
  551D2304 18301680 14AAD863 8DE1B749 928AFDD8 7A3FC27B C1838544 38301D06 
  03551D0E 04160414 AAD8638D E1B74992 8AFDD87A 3FC27BC1 83854438 300D0609 
  2A864886 F70D0101 05050003 81810003 2AC29609 C3C90680 59AD6453 FAF00D55 
  C7A1854C A19D5EFE FF023634 D0CFE672 D497425F 445F297E 5F0EB106 060EB438 
  DD7627F1 855EAC71 46F07426 9D5A4171 5F6E5853 E1C05531 3C32CFFB 3DB3455D 
  83034B46 F1C3B752 40365A3C 133BD49D 1E05858F A393CFD4 5EBE1CCE D5E1D1D3 
  DDA0A8A3 23311A87 905B88B9 56DD91
  	quit
!
!
!
!
!
!
spanning-tree mode mst
spanning-tree extend system-id
!
spanning-tree mst configuration
 name REGION-PELLET
 revision 1
 instance 1 vlan 600, 650, 700
 instance 2 vlan 120, 130, 140, 150, 160
!
!
vlan internal allocation policy ascending
!         
ip ssh time-out 60
ip ssh version 2
!
!
!
!
!
interface Port-channel1
 description Vers_SWR1_Principal
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
!
interface Port-channel2
 description Vers_SWR2_Secours
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
!
interface FastEthernet0/1
 description TEST_COMPTA_V120
 switchport access vlan 120
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/2
 description TEST_COMPTA_V120
 switchport access vlan 120
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/3
 description TEST_COMPTA_V120
 switchport access vlan 120
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/4
 description TEST_COMPTA_V120
 switchport access vlan 120
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/5
 description TEST_TERRAIN_V130
 switchport access vlan 130
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/6
 description TEST_TERRAIN_V130
 switchport access vlan 130
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/7
 description TEST_TERRAIN_V130
 switchport access vlan 130
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/8
 description TEST_TERRAIN_V130
 switchport access vlan 130
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/9
 description TEST_ADMIN_V140
 switchport access vlan 140
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/10
 description TEST_ADMIN_V140
 switchport access vlan 140
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/11
 description TEST_ADMIN_V140
 switchport access vlan 140
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/12
 description TEST_ADMIN_V140
 switchport access vlan 140
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/13
 description TEST_PROD_V150
 switchport access vlan 150
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/14
 description TEST_PROD_V150
 switchport access vlan 150
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/15
 description TEST_PROD_V150
 switchport access vlan 150
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/16
 description TEST_PROD_V150
 switchport access vlan 150
 switchport mode access
 spanning-tree portfast
!         
interface FastEthernet0/17
 description TEST_LOGIS_V160
 switchport access vlan 160
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/18
 description TEST_LOGIS_V160
 switchport access vlan 160
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/19
 description TEST_ADMIN_INFO_V600
 switchport access vlan 600
 switchport mode access
 spanning-tree portfast
!
interface FastEthernet0/20
!
interface FastEthernet0/21
 description Uplink_SWR2
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
 channel-group 2 mode active
!
interface FastEthernet0/22
 description Uplink_SWR2
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
 channel-group 2 mode active
!
interface FastEthernet0/23
 description Uplink_SWR1
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
 channel-protocol lacp
 channel-group 1 mode active
!
interface FastEthernet0/24
 description Uplink_SWR1
 switchport trunk native vlan 999
 switchport trunk allowed vlan 120-160,600,650,700,999
 switchport mode trunk
 channel-protocol lacp
 channel-group 1 mode active
!
interface GigabitEthernet0/1
!
interface GigabitEthernet0/2
!
interface Vlan1
 no ip address
 shutdown
!
interface Vlan600
 description Management_SW2
 ip address 172.16.0.11 255.255.255.0
!
ip default-gateway 172.16.0.254
ip http server
ip http secure-server
access-list 10 permit 172.16.0.0 0.0.0.255
!
vstack
!         
line con 0
line vty 0 4
 access-class 10 in
 exec-timeout 5 0
 login local
 transport input ssh
line vty 5 15
 access-class 10 in
 exec-timeout 5 0
 login local
 transport input ssh
!
end
