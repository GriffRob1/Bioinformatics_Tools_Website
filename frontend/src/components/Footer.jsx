import React from 'react';

export default function Footer() {
    return (
        <footer className={'container footer'}>
            <div className={'copyright-information'}>Copyright 2025 by Big Griff LLC. All Rights Reserved.</div>
            <div className={'container social-media'}>
                <p>Follow us on social media!</p>
                <img src={'/images/317727_facebook_social media_social_icon.svg'}  alt={'facebook'}/>
                <img src={'/images/1298747_instagram_brand_logo_social media_icon.svg'}  alt={'instagram'}/>
                <img src={'/images/4202098_github_code_developer_logo_icon.svg'}  alt={'github'}/>
            </div>
        </footer>
    )
}