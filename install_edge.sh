#! /bin/bash

token="e3a9c925b78599907fa0ac5b10d3e5bc41644226bf7cbda2df50f0a5954cae90.eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJleHAiOjE3MTU1ODgzMDh9.jJF4LdUcel7zp-7ILXgHTVXFwVFNu6T42ScXeC9uRXw"

# Pre-requisites for installing container runtime
echo "Setting Kernel parameter for ipv4 packet forwarding"
cat <<EOF | tee /etc/sysctl.d/k8s.conf
net.ipv4.ip_forward = 0
EOF
sysctl --system

# Installing container runtime interface plugin - containerd
echo "Installing containerd as the container runtime"
apt-get install containerd
mkdir -p /etc/containerd
cp config.toml /etc/containerd/config.toml
systemctl restart containerd


wget https://github.com/containernetworking/plugins/releases/download/v1.3.0/cni-plugins-linux-amd64-v1.3.0.tgz
sudo mkdir -p /opt/cni/bin
sudo tar Cxzvf /opt/cni/bin cni-plugins-linux-amd64-v1.3.0.tgz
rm cni-plugins-*
sudo cat << EOF | sudo tee /etc/cni/net.d/10-containerd-net.conflist
{
  "cniVersion": "1.0.0",
  "name": "containerd-net",
  "plugins": [
    {
      "type": "bridge",
      "bridge": "cni0",
      "isGateway": true,
      "ipMasq": true,
      "promiscMode": true,
      "ipam": {
        "type": "host-local",
        "ranges": [
          [{
            "subnet": "10.88.0.0/16"
          }],
          [{
            "subnet": "2001:4860:4860::/64"
          }]
        ],
        "routes": [
          { "dst": "0.0.0.0/0" },
          { "dst": "::/0" }
        ]
      }
    },
    {
      "type": "portmap",
      "capabilities": {"portMappings": true}
    }
  ]
}
EOF
sudo systemctl restart containerd


# Installing edgecore service
echo "Installing edgecore service"
cd /home/aadesh/install_edge
sed -i -e "s|token: .*|token: ${token}|g" edgecore.yaml
mkdir -p /etc/kubeedge/config
cp edgecore.yaml /etc/kubeedge/config
cd /opt
wget https://github.com/kubeedge/kubeedge/releases/download/v1.14.5/kubeedge-v1.14.5-linux-amd64.tar.gz
tar -zxvf kubeedge-v1.14.5-linux-amd64.tar.gz
cp kubeedge-v1.14.5-linux-amd64/edge/edgecore /usr/local/bin/edgecore
cat <<EOF > /etc/systemd/system/edgecore.service
[Unit]
Description=edgecore.service

[Service]
Type=Simple
ExecStart=/usr/local/bin/edgecore
Restart=always
RestartSec=10
Environment=DEPLOY_MQTT_CONTAINER=TRUE
KillMode=process

[Install]
WantedBy=multi-user.target
EOF
systemctl daemon-reload
systemctl enable --now edgecore
systemctl start edgecore
